use crate::traits::InputMF;
use std::cmp::Ordering;

#[derive(Debug, PartialEq)]
struct TriMF {
    slope1: f32,
    intcpt1: f32,
    slope2: f32,
    intcpt2: f32,
    a: f32,
    xstar: f32,
    b: f32,
}

impl TriMF {
    fn new(a: f32, xstar: f32, b: f32) -> Self {
        let mut a = a;
        let mut xstar = xstar;
        let mut b = b;
        if (a > xstar) || (xstar > b) || (b < a) {
            let mut params = [a, xstar, b];
            params.sort_unstable_by(f32::total_cmp);
            [a, xstar, b] = params;
        }

        let slope1 = if (xstar - a).abs() > f32::EPSILON {
            1.0 / (xstar - a)
        } else {
            0.0
        };
        let slope2 = if (b - xstar).abs() > f32::EPSILON {
            -1.0 / (b - xstar)
        } else {
            0.0
        };

        let intcpt1 = -slope1 * a;
        let intcpt2 = -slope2 * b;
        Self {
            slope1,
            intcpt1,
            slope2,
            intcpt2,
            a,
            xstar,
            b,
        }
    }
}

impl InputMF for TriMF {
    fn fuzzify(&self, x: f32) -> f32 {
        let x = x.clamp(self.a, self.b);
        match x.total_cmp(&self.xstar) {
            Ordering::Less => x.mul_add(self.slope1, self.intcpt1).clamp(0.0, 1.0),
            Ordering::Greater => x.mul_add(self.slope2, self.intcpt2).clamp(0.0, 1.0),
            Ordering::Equal => 1.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::TriMF;
    use crate::traits::InputMF;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_trimf_new() {
        let isosc012 = TriMF {
                slope1: 1.0,
                intcpt1: 0.0,
                slope2: -1.0,
                intcpt2: 2.0,
                a: 0.0,
                xstar: 1.0,
                b: 2.0,
            };
        let desc_ramp = TriMF {
                slope1: 0.0,
                intcpt1: 0.0,
                slope2: -1.0,
                intcpt2: 1.0,
                a: 0.0,
                xstar: 0.0,
                b: 1.0,
            };
        let asc_ramp = TriMF {
                slope1: 1.0,
                intcpt1: 0.0,
                slope2: 0.0,
                intcpt2: 0.0,
                a: 0.0,
                xstar: 1.0,
                b: 1.0,
            };

        // Test well-formed parameter input
        assert_eq!(isosc012,
            TriMF::new(0.0, 1.0, 2.0)
        );
        assert_eq!(
            desc_ramp,
            TriMF::new(0.0, 0.0, 1.0)
        );
        assert_eq!(
            asc_ramp,
            TriMF::new(0.0, 1.0, 1.0)
        );

        // Test unsorted params
        assert_eq!(isosc012,
            TriMF::new(1.0, 2.0, 0.0)
        );
        assert_eq!(isosc012,
            TriMF::new(1.0, 0.0, 2.0)
        );
        assert_eq!(
            desc_ramp,
            TriMF::new(0.0, 1.0, 0.0)
        );
        assert_eq!(
            desc_ramp,
            TriMF::new(1.0, 0.0, 0.0)
        );
        assert_eq!(
            asc_ramp,
            TriMF::new(1.0, 1.0, 0.0)
        );
        assert_eq!(
            asc_ramp,
            TriMF::new(1.0, 0.0, 1.0)
        );
    }

    #[test]
    fn test_trimf_fuzzify() {
        let isosc012 = TriMF::new(0.0, 1.0, 2.0);

        assert_approx_eq!(isosc012.fuzzify(-1.0), 0.0);
        assert_approx_eq!(isosc012.fuzzify(0.0), 0.0);
        assert_approx_eq!(isosc012.fuzzify(0.25), 0.25);
        assert_approx_eq!(isosc012.fuzzify(0.5), 0.5);
        assert_approx_eq!(isosc012.fuzzify(0.75), 0.75);
        assert_approx_eq!(isosc012.fuzzify(1.0), 1.0);
        assert_approx_eq!(isosc012.fuzzify(1.25), 0.75);
        assert_approx_eq!(isosc012.fuzzify(1.5), 0.5);
        assert_approx_eq!(isosc012.fuzzify(1.75), 0.25);
        assert_approx_eq!(isosc012.fuzzify(2.0), 0.0);
        assert_approx_eq!(isosc012.fuzzify(3.0), 0.0);
    }

    #[test]
    fn test_descramp_mf_fuzzify() {
        let desc_ramp = TriMF::new(0.0, 0.0, 1.0);

        assert_approx_eq!(desc_ramp.fuzzify(-1.0), 1.0);
        assert_approx_eq!(desc_ramp.fuzzify(0.0), 1.0);
        assert_approx_eq!(desc_ramp.fuzzify(0.25), 0.75);
        assert_approx_eq!(desc_ramp.fuzzify(0.5), 0.5);
        assert_approx_eq!(desc_ramp.fuzzify(0.75), 0.25);
        assert_approx_eq!(desc_ramp.fuzzify(1.0), 0.0);
        assert_approx_eq!(desc_ramp.fuzzify(2.0), 0.0);
    }

    #[test]
    fn test_ascramp_mf_fuzzify() {
        let asc_ramp = TriMF::new(0.0, 1.0, 1.0);

        assert_approx_eq!(asc_ramp.fuzzify(-1.0), 0.0);
        assert_approx_eq!(asc_ramp.fuzzify(0.0), 0.0);
        assert_approx_eq!(asc_ramp.fuzzify(0.25), 0.25);
        assert_approx_eq!(asc_ramp.fuzzify(0.5), 0.5);
        assert_approx_eq!(asc_ramp.fuzzify(0.75), 0.75);
        assert_approx_eq!(asc_ramp.fuzzify(1.0), 1.0);
        assert_approx_eq!(asc_ramp.fuzzify(2.0), 1.0);
    }

    #[test]
    fn test_floats_fuzzify() {
        let asc_ramp = TriMF::new(0.0, 100.0, 100.0);

        assert_approx_eq!(asc_ramp.fuzzify(-f32::INFINITY), 0.0);
        assert_approx_eq!(asc_ramp.fuzzify(f32::EPSILON.powi(2)), 0.0);
        assert_approx_eq!(asc_ramp.fuzzify(f32::EPSILON), 0.0);
        assert_approx_eq!(asc_ramp.fuzzify(f32::INFINITY), 1.0);
        assert!(asc_ramp.fuzzify(-f32::NAN).is_nan());
        assert!(asc_ramp.fuzzify(f32::NAN).is_nan());
        assert_approx_eq!(asc_ramp.fuzzify(-0.0), 0.0);
        assert_approx_eq!(asc_ramp.fuzzify(1.0e-5), 1.0e-7, f32::EPSILON);

        let mf = TriMF::new(0.0, 6.0e8, 8.0e8);
        assert_approx_eq!(mf.fuzzify(7.0e8 + 3.0e-3), 0.5 + 3.0e-12, f32::EPSILON)
    }
}
