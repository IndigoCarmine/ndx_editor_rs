use crate::error::{NdxError, Result};
use crate::expr::{EvalCtx, auto_name, parse_expr};
use crate::model::{Group, IndexFile};
use crate::system::SystemCtx;
use crate::universe::UniverseSpec;

pub struct Selection {
    /// Index of the group that was appended.
    pub id: usize,
    pub name: String,
    pub len: usize,
    pub warnings: Vec<String>,
}

/// Evaluate `src` against `ndx` and append the result as a new group.
///
/// An empty result is a warning, not an error, unless `error_on_empty`: `A & B` coming out empty
/// is frequently the answer the user was looking for, and make_ndx creates the group too.
pub fn select(
    ndx: &mut IndexFile,
    src: &str,
    name: Option<&str>,
    spec: &UniverseSpec,
    system: &SystemCtx,
    error_on_empty: bool,
) -> Result<Selection> {
    let expr = parse_expr(src)?;

    let (atoms, mut warnings) = {
        let cx = EvalCtx::new(ndx, spec, system);
        let set = cx.run(&expr)?;
        (set.as_slice().to_vec(), cx.take_warnings())
    };

    if atoms.is_empty() {
        if error_on_empty {
            return Err(NdxError::EmptyResult);
        }
        warnings.push(format!("the expression {src:?} selected no atoms"));
    }

    let name = match name {
        Some(n) => n.to_string(),
        None => ndx.unique_name(&auto_name(ndx, &expr)),
    };
    let len = atoms.len();
    let id = ndx.push(Group::new(name.clone(), atoms));

    Ok(Selection {
        id,
        name,
        len,
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Group;

    fn ndx() -> IndexFile {
        IndexFile {
            groups: vec![
                Group::new("System", (1..=10).collect()),
                Group::new("Protein", vec![1, 2, 3]),
                Group::new("SOL", vec![3, 4, 5]),
            ],
        }
    }

    fn sel(f: &mut IndexFile, src: &str) -> Result<Selection> {
        select(
            f,
            src,
            None,
            &UniverseSpec::Auto,
            &SystemCtx::default(),
            false,
        )
    }

    #[test]
    fn appends_a_named_group() {
        let mut f = ndx();
        let s = sel(&mut f, "1 & !2").unwrap();
        assert_eq!(s.id, 3);
        assert_eq!(s.name, "Protein_&_!SOL");
        assert_eq!(f.groups[3].atoms, [1, 2]);
    }

    #[test]
    fn explicit_name_wins() {
        let mut f = ndx();
        let s = select(
            &mut f,
            "1 | 2",
            Some("Both"),
            &UniverseSpec::Auto,
            &SystemCtx::default(),
            false,
        )
        .unwrap();
        assert_eq!(s.name, "Both");
    }

    #[test]
    fn empty_result_warns_but_creates() {
        let mut f = ndx();
        let s = sel(&mut f, "1 & !1").unwrap();
        assert_eq!(s.len, 0);
        assert_eq!(f.len(), 4);
        assert!(s.warnings.iter().any(|w| w.contains("no atoms")));
    }

    #[test]
    fn error_on_empty_opts_out() {
        let mut f = ndx();
        let r = select(
            &mut f,
            "1 & !1",
            None,
            &UniverseSpec::Auto,
            &SystemCtx::default(),
            true,
        );
        assert!(matches!(r, Err(NdxError::EmptyResult)));
        assert_eq!(f.len(), 3, "nothing should have been appended");
    }

    #[test]
    fn a_repeated_auto_name_gets_a_suffix() {
        let mut f = ndx();
        sel(&mut f, "1 & 2").unwrap();
        let s = sel(&mut f, "1 & 2").unwrap();
        assert_eq!(s.name, "Protein_&_SOL_2");
    }
}
