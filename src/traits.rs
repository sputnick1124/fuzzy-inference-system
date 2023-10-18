use std::ops::Add;

pub trait Membership {}

pub trait InputMF {
    fn fuzzify(&self, input: f32) -> f32;
}

pub trait OutputMF {
    type Output: Membership;

    fn membership(&self, activation: f32) -> Self::Output;
}
