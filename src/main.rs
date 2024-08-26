//! Main file exampling setup and parameters
#![allow(warnings)]
#![allow(non_snake_case)]
#![allow(warnings)]
#![feature(iter_next_chunk)]
use rust_elgamal::{Scalar, Ciphertext};
use rand_chacha::ChaCha20Rng;
use rand::SeedableRng;

use std::time::Duration;

use merlin::Transcript;

use crate::utils::enums::EGInp;

mod arguers;
mod traits;
mod utils;
mod provers;

use crate::traits::traits::{HeapSize, Timeable};
use crate::utils::*;
use crate::provers::*;

fn run_shuffle_argument(m: usize, n: usize, mu: usize) -> (u128, u128, usize) {

    let N = m * n;

    assert!(m % mu == 0);

    //Setup rng and common reference key(public key for ElGamal,
    // commitment key for Pedersen)
    let mut rng = ChaCha20Rng::from_entropy();
    let mut cr = arguers::CommonRef::new(N as u64, rng);

    //Create Open Deck from scalars
    let deck: Vec<Scalar> = (0..N).map(|card| Scalar::from(card as u64))
                                    .collect();
    //Get blinding values for the deck
    let deck_r: Vec<Scalar> = (0..N).map(|_| cr.rand_scalar()).collect();
    //Encrypt the decks in ElGamal
    let C_deck: Vec<Ciphertext> = deck.iter()
                                        .zip(deck_r.clone())
                                        .map(|(card, r)| 
                                             cr.encrypt(&EGInp::Scal(card.clone()), &r)
                                        )
                                .collect();

    //Get a random permutation
    let permutation: Vec<u64> = cr.rand_perm(&(1..=(N as u64)).collect());

    //Get more blinding values for hiding the permuted cards
    let rho: Vec<Scalar> = (0..N).map(|_| cr.rand_scalar()).collect();

    //Permute the cards while hiding them with another ElGamal Encryption
    let C_permd: Vec<Ciphertext> = permutation.iter()
                                                .zip(rho.iter())
                                                .map(|(pi, r)| {
                                                    let card: Ciphertext = C_deck[*pi as usize - 1].clone();
                                                    card + cr.encrypt(&EGInp::Scal(Scalar::zero()), &r)

                                                    }
                                                )
                                                .collect();

    //Create the main shuffle prover
    let mut prover_transcript = Transcript::new(b"ShuffleProof");
    let mut shuffle_prover = prover::ShuffleProver::new(
                            m,
                            n,
                            mu,
                            C_deck,
                            C_permd,
                            permutation,
                            rho,
                            cr
                        );
                                
    let prove_time = shuffle_prover.start_time();
    //Recieve the provers used in the argument and the final proof
    let (mut zero_prover,
         mut sv_prover,
         mut hadam_prover,
         mut prod_prover,
         mut mexp_prover,
         mut shuffle_proof) = shuffle_prover.prove(&mut prover_transcript);

    let prove_time = shuffle_prover.elapsed(prove_time);
    let proof_size = shuffle_proof.heap_size();

    //Create a clean transcript for the verification
    let mut verifier_transcript = Transcript::new(b"ShuffleProof");

    let verify_time = shuffle_prover.start_time();
    //Assert verification returns Ok
    assert!(shuffle_prover
            .verify(&mut verifier_transcript, shuffle_proof, 
                    zero_prover,
                    sv_prover,
                    hadam_prover,
                    prod_prover,
                    mexp_prover)
            .is_ok());
    let verify_time = shuffle_prover.elapsed(verify_time);

    (prove_time, verify_time, proof_size)
}


fn main() -> Result<(), Box<dyn std::error::Error>>{
    use plotters::prelude::*;
    use std::collections::HashMap;
    use itertools::iproduct;
    use signal_hook::flag;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, Ordering};

    let term = Arc::new(AtomicBool::new(false));
    flag::register(signal_hook::consts::SIGINT, Arc::clone(&term))?;

    let factor = 1;
    let mu: usize = 1;

    let trials = 3;

    let m: Vec<usize> = (2..trials).collect();
    let n: Vec<usize> = m.iter().map(|a| factor * a).collect();


    let mut max_time: u128 = 0;
    let mut max_size: usize = 0;
    let mut indexes: Vec<usize> = vec![];
    let mut p_time: Vec<u128> = vec![];
    let mut v_time: Vec<u128> = vec![];
    let mut p_size: Vec<usize> = vec![];

    let mut it = m.into_iter().zip(n.into_iter());
    let mut round = 0;
    while !term.load(Ordering::Relaxed) {
        let (m_, n_) = it.next().unwrap();
        println!("{}", m_ * n_);
        let metric = run_shuffle_argument(m_, n_, mu);
        if metric.0 > max_time {max_time = metric.0;}
        if metric.1 > max_time {max_time = metric.1;}
        if metric.2 > max_size {max_size = metric.2;}
        indexes.push(m_ * n_);
        p_time.push(metric.0);
        v_time.push(metric.1);
        p_size.push(metric.2);

        round += 1;
        if round >= trials {break;}
    }

    println!("DRAWING");

    let root = BitMapBackend::new("time.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Time metrics in ms / Card Count", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(50)
        .build_cartesian_2d(0..(*indexes.last().unwrap()), 0..max_time)?;

    chart.configure_mesh().draw()?;

    chart.draw_series(LineSeries::new(
            indexes.clone().into_iter().zip(p_time.into_iter()),
            Into::<ShapeStyle>::into(&RED).stroke_width(2)
            ))?
        .label("proof time")
        .legend(|(x, y)| PathElement::new(vec![(x,y), (x+20, y)], &RED));
    chart.draw_series(LineSeries::new(
            indexes.clone().into_iter().zip(v_time.into_iter()),
            Into::<ShapeStyle>::into(&BLUE).stroke_width(2)
            ))?
        .label("verify time")
        .legend(|(x, y)| PathElement::new(vec![(x,y), (x+20, y)], &BLUE));

    chart
        .configure_series_labels()
        .draw()?;

    root.present()?;

    let root = BitMapBackend::new("space.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Proof Size in bytes / Card Count", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(50)
        .build_cartesian_2d(0..(*indexes.last().unwrap()), 0..max_size)?;

    chart.configure_mesh().draw()?;

    chart.draw_series(LineSeries::new(
            indexes.into_iter().zip(p_size.into_iter()),
            Into::<ShapeStyle>::into(&GREEN).stroke_width(2)
            ))?
        .label("proof size")
        .legend(|(x, y)| PathElement::new(vec![(x,y), (x+20, y)], &GREEN));
    chart
        .configure_series_labels()
        .draw()?;

    root.present()?;

    Ok(())
}
