extern crate threshold_secret_sharing as tss;

fn main() {

    let ref pss = tss::PSS_4_26_3;
    println!("\
    Using parameters that: \n \
     - allow {} values to be packed together \n \
     - give a security threshold of {} \n \
     - require {} of the {} shares to reconstruct in the basic case",
        pss.secret_count,
        pss.threshold,
        pss.reconstruct_limit,
        pss.share_count
    );

    // define inputs
    let secrets_1 = vec![1, 2, 3];
    println!("\nFirst input vector:  {:?}", &secrets_1);
    let secrets_2 = vec![4, 5, 6];
    println!("Second input vector: {:?}", &secrets_2);

    // secret share inputs
    let shares_1 = pss.share(&secrets_1);
    println!("\nSharing of first vector gives random shares:\n{:?}", &shares_1);
    let shares_2 = pss.share(&secrets_2);
    println!("\nSharing of second vector gives random shares:\n{:?}", &shares_2);

    // in the following, 'positivise' is used to map (potentially negative)
    // values to their equivalent positive representation in Z_p for usability
    use tss::numtheory::positivise;

    // add shares point-wise
    let sum_shares: Vec<i64> = shares_1.iter().zip(&shares_2).map(|(a, b)| (a + b) % pss.prime).collect();
    println!("\nAdding the shares point-wise gives shares:\n{:?}", &sum_shares);

    // reconstruct sum, using same reconstruction limit
    let sum_reconstruct_limit = pss.reconstruct_limit;
    let sum_indices: Vec<usize> = (0..sum_reconstruct_limit).collect();
    let sum_shares = &sum_shares[0..sum_reconstruct_limit];
    let sum_secrets = pss.reconstruct(&sum_indices, sum_shares);
    println!(
        "... and reconstructing using {} of these gives output vector: {:?}",
        sum_reconstruct_limit,
        positivise(&sum_secrets, pss.prime)
    );

    // multiply shares point-wise
    let product_shares: Vec<i64> = shares_1.iter().zip(&shares_2).map(|(a, b)| (a * b) % pss.prime).collect();
    println!("\nMultiplying the shares point-wise gives shares:\n{:?}", &product_shares);

    // reconstruct product, using double reconstruction limit (minus one)
    let product_reconstruct_limit = pss.reconstruct_limit*2-1;
    let product_indices: Vec<usize> = (0..product_reconstruct_limit).collect();
    let product_shares = &product_shares[0..product_reconstruct_limit];
    let product_secrets = pss.reconstruct(&product_indices, product_shares);
    println!(
        "... and reconstructing using {} of these gives output vector: {:?}",
        product_reconstruct_limit,
        positivise(&product_secrets, pss.prime)
    );

}
