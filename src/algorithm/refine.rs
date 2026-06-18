//! Algorithm 2 - Iterative Refinement: solve eq. 40 on `T^est`, drop slack
//! times, add violated local maxima, until convergence.

/// Values within `PLATEAU_EPS` of each other are treated as a flat top.
const PLATEAU_EPS: f64 = 1e-12;

/// Indices of `g` that are local maxima **and** exceed `threshold`.
///
/// A flat top (run of values within [`PLATEAU_EPS`]) yields a single
/// representative (the plateau midpoint). Endpoints are local maxima by the
/// boundary rule (compared only against their one in-bounds neighbour). The
/// global maximum is always a local maximum, so a violated global max is always
/// returned — guaranteeing Algorithm 2 makes progress.
#[allow(dead_code)] // wired into refine() in Task 3
pub(super) fn violated_local_maxima(g: &[f64], threshold: f64) -> Vec<usize> {
    let n = g.len();
    let mut out = Vec::new();
    let mut k = 0usize;
    while k < n {
        // Extent of the flat run [k..=j] of values ~equal to g[k].
        let mut j = k;
        while j + 1 < n && (g[j + 1] - g[k]).abs() <= PLATEAU_EPS {
            j += 1;
        }
        let left_ok = k == 0 || g[k - 1] < g[k];
        let right_ok = j == n - 1 || g[j + 1] < g[k];
        if left_ok && right_ok && g[k] > threshold {
            out.push((k + j) / 2); // plateau midpoint representative
        }
        k = j + 1;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interior_peak() {
        assert_eq!(violated_local_maxima(&[0.0, 2.0, 0.0], 1.0), vec![1]);
    }

    #[test]
    fn left_endpoint_peak() {
        assert_eq!(violated_local_maxima(&[3.0, 1.0, 0.5], 1.0), vec![0]);
    }

    #[test]
    fn right_endpoint_peak() {
        assert_eq!(violated_local_maxima(&[0.5, 1.0, 3.0], 1.0), vec![2]);
    }

    #[test]
    fn flat_top_two_yields_one_midpoint() {
        // plateau [1,2] -> midpoint (1+2)/2 = 1.
        assert_eq!(violated_local_maxima(&[0.0, 2.0, 2.0, 0.0], 1.0), vec![1]);
    }

    #[test]
    fn flat_top_three_yields_one_midpoint() {
        // plateau [1,3] -> midpoint (1+3)/2 = 2.
        assert_eq!(
            violated_local_maxima(&[0.0, 2.0, 2.0, 2.0, 0.0], 1.0),
            vec![2]
        );
    }

    #[test]
    fn monotone_increasing_picks_last() {
        assert_eq!(violated_local_maxima(&[0.0, 1.0, 2.0, 3.0], 1.0), vec![3]);
    }

    #[test]
    fn monotone_decreasing_picks_first() {
        assert_eq!(violated_local_maxima(&[3.0, 2.0, 1.0, 0.0], 1.0), vec![0]);
    }

    #[test]
    fn all_below_threshold_is_empty() {
        assert!(violated_local_maxima(&[0.1, 0.2, 0.1], 1.0).is_empty());
    }

    #[test]
    fn two_separated_peaks() {
        assert_eq!(
            violated_local_maxima(&[0.0, 2.0, 0.5, 3.0, 0.0], 1.0),
            vec![1, 3]
        );
    }

    #[test]
    fn threshold_filters_low_peak() {
        // peak at idx 1 (1.5 > 1) kept; peak at idx 3 (0.8 <= 1) dropped.
        assert_eq!(
            violated_local_maxima(&[0.0, 1.5, 0.0, 0.8, 0.0], 1.0),
            vec![1]
        );
    }

    #[test]
    fn single_element_above_and_below() {
        assert_eq!(violated_local_maxima(&[5.0], 1.0), vec![0]);
        assert!(violated_local_maxima(&[0.5], 1.0).is_empty());
    }
}
