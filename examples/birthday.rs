use itertools::Itertools;
use phantom_zone::*;
use rand::{thread_rng, RngCore};
use rayon::prelude::*;
use std::array;

const YEAR_ZERO: usize = 1900;
const YEAR_MAX: usize = 2040;
const DAYS_MONTH: [usize; 12] = [31, 30, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

const DATE_ID_BITS: usize = 16;
const DATA_BITS: usize = 64;

// 1900-01-01: date 0
// 2040-01-01: date 54900
// 16 bit date
fn date_to_id(year: usize, month: usize, day: usize) -> u64 {
    // Input validation
    assert!(YEAR_ZERO <= year && year < YEAR_MAX);
    let year_rel = year - YEAR_ZERO;
    assert!(1 <= month && month <= 12);
    let month_rel = month - 1;
    assert!(1 <= day && day <= DAYS_MONTH[month - 1]);
    let day_rel = day - 1;

    let id = year_rel * 366 + (0..month_rel).map(|m| DAYS_MONTH[m]).sum::<usize>() + day_rel;
    assert!(id < 2usize.pow(DATE_ID_BITS as u32)); // id should be DATE_ID_BITS bits
    id as u64
}

fn u64_to_binary<const N: usize>(v: u64) -> [bool; N] {
    assert!((v as u128) < 2u128.pow(N as u32));
    let mut result = [false; N];
    for i in 0..N {
        if (v >> i) & 1 == 1 {
            result[i] = true;
        }
    }
    result
}

fn binary_to_u64(v: &[bool]) -> u64 {
    assert!(v.len() <= 64);
    let mut result = 0;
    for i in 0..v.len() {
        if v[i] {
            result |= 1 << i;
        }
    }
    result
}

fn birthday_match(
    date_id_a: &[bool; DATE_ID_BITS],
    data_a: &[bool; DATA_BITS],
    date_id_b: &[bool; DATE_ID_BITS],
    data_b: &[bool; DATA_BITS],
) -> ([bool; DATA_BITS], [bool; DATA_BITS]) {
    let not_is_match = date_id_a
        .iter()
        .zip(date_id_b.iter())
        .map(|(bit_a, bit_b)| bit_a ^ bit_b)
        .reduce(|acc, r| &acc | &r)
        .unwrap();
    let is_match = !&not_is_match;
    let masked_data_a = array::from_fn(|i| &data_a[i] & &is_match);
    let masked_data_b = array::from_fn(|i| &data_b[i] & &is_match);
    (masked_data_a, masked_data_b)
}

fn birthday_match_fhe(
    date_id_a: &[FheBool; DATE_ID_BITS],
    data_a: &[FheBool; DATA_BITS],
    date_id_b: &[FheBool; DATE_ID_BITS],
    data_b: &[FheBool; DATA_BITS],
) -> ([FheBool; DATA_BITS], [FheBool; DATA_BITS]) {
    let not_is_match = date_id_a
        .iter()
        .zip(date_id_b.iter())
        .map(|(bit_a, bit_b)| bit_a ^ bit_b)
        .reduce(|acc, r| &acc | &r)
        .unwrap();
    let is_match = !&not_is_match;
    let masked_data_a = array::from_fn(|i| &data_a[i] & &is_match);
    let masked_data_b = array::from_fn(|i| &data_b[i] & &is_match);
    (masked_data_a, masked_data_b)
}

fn birthday_match_fhe_par(
    date_id_a: &[FheBool; DATE_ID_BITS],
    data_a: &[FheBool; DATA_BITS],
    date_id_b: &[FheBool; DATE_ID_BITS],
    data_b: &[FheBool; DATA_BITS],
) -> ([FheBool; DATA_BITS], [FheBool; DATA_BITS]) {
    rayon::ThreadPoolBuilder::new()
        .build_scoped(
            // Initialize thread-local storage parameters
            |thread| {
                set_parameter_set(ParameterSelector::NonInteractiveLTE2Party);
                thread.run()
            },
            // Run parallel code under this pool
            |pool| pool.install(|| _birthday_match_fhe_par(date_id_a, data_a, date_id_b, data_b)),
        )
        .unwrap()
}

fn _birthday_match_fhe_par(
    date_id_a: &[FheBool; DATE_ID_BITS],
    data_a: &[FheBool; DATA_BITS],
    date_id_b: &[FheBool; DATE_ID_BITS],
    data_b: &[FheBool; DATA_BITS],
) -> ([FheBool; DATA_BITS], [FheBool; DATA_BITS]) {
    let not_is_match = date_id_a
        .par_iter()
        .zip(date_id_b.par_iter())
        .map(|(bit_a, bit_b)| bit_a ^ bit_b)
        .reduce_with(|acc, r| &acc | &r)
        .unwrap();
    let is_match = !&not_is_match;
    let masked_data_a: Vec<_> = data_a.par_iter().map(|byte| byte & &is_match).collect();
    let masked_data_b: Vec<_> = data_b.par_iter().map(|byte| byte & &is_match).collect();
    (
        masked_data_a.try_into().unwrap_or_else(|_| panic!()),
        masked_data_b.try_into().unwrap_or_else(|_| panic!()),
    )
}

fn main() {
    set_parameter_set(ParameterSelector::NonInteractiveLTE2Party);

    // set application's common reference seed
    let mut seed = [0u8; 32];
    thread_rng().fill_bytes(&mut seed);
    set_common_reference_seed(seed);

    // let n_users = 100;
    let no_of_parties = 2;

    // Clide side //

    // Generate client keys
    let cks = (0..no_of_parties).map(|_| gen_client_key()).collect_vec();
    let ck_a = &cks[0];
    let ck_b = &cks[1];

    let date_a = (1989, 03, 25);
    let data_a = 0x1234;
    let date_b = (1989, 03, 26);
    let data_b = 0x5678;

    let date_id_a = date_to_id(date_a.0, date_a.1, date_a.2);
    let date_id_b = date_to_id(date_b.0, date_b.1, date_b.2);
    let date_id_a_bin: [bool; DATE_ID_BITS] = u64_to_binary(date_id_a);
    let date_id_b_bin: [bool; DATE_ID_BITS] = u64_to_binary(date_id_b);
    let data_a_bin: [bool; DATA_BITS] = u64_to_binary(data_a);
    let data_b_bin: [bool; DATA_BITS] = u64_to_binary(data_b);

    // Native example
    // let res = birthday_match(&date_id_a_bin, &data_a_bin, &date_id_b_bin,
    // &data_b_bin); println!("{:?}", res);

    // client A encrypts its private inputs
    //
    // Clients encrypt their private inputs in a seeded batched ciphertext using
    // their private RLWE secret `u_j`.
    let input_a: Vec<bool> = date_id_a_bin
        .iter()
        .chain(data_a_bin.iter())
        .copied()
        .collect();

    // client B encrypts its private inputs
    //
    // Clients encrypt their private inputs in a seeded batched ciphertext using
    // their private RLWE secret `u_j`.
    let input_b: Vec<bool> = date_id_b_bin
        .iter()
        .chain(data_b_bin.iter())
        .copied()
        .collect();

    // Clients independently generate their server key shares
    //
    // We assign user_id 0 to client A, user_id 1 to client B.
    //
    // Note that `user_id`s must be unique among the clients and must be less than
    // total number of clients.
    let now = std::time::Instant::now();
    let server_key_shares = cks
        .iter()
        .enumerate()
        .map(|(id, k)| gen_server_key_share(id, no_of_parties, k))
        .collect_vec();
    println!("Clients server key share gen time: {:?}", now.elapsed());

    // Each client uploads their server key shares and encrypted private inputs to
    // the server in a single shot message.

    // Server side //

    // Server receives server key shares from each client and proceeds to aggregate
    // them to produce the server key. After this point, server can use the server
    // key to evaluate any arbitrary function on encrypted private inputs from
    // the fixed set of clients

    // aggregate server shares and generate the server key
    let now = std::time::Instant::now();
    let server_key = aggregate_server_key_shares(&server_key_shares);
    server_key.set_server_key();
    println!("Server key gen time: {:?}", now.elapsed());

    // Server proceeds to extract private inputs sent by clients
    //
    // To extract client A's (with user_id=0) private inputs we first key switch
    // client A's private inputs from theit secret `u_j` to ideal secret of the
    // mpc protocol. To indicate we're key switching client A's private
    // input we supply client A's `user_id` i.e. we call `key_switch(0)`.
    // Then we extract the first ciphertext by calling `extract_at(0)`.
    //
    // Since client A only encrypts 1 input in batched ciphertext, calling
    // extract_at(index) for `index` > 0 will panic. If client A had more
    // private inputs then we can either extract them all at once with
    // `extract_all` or first `many` of them with `extract_many(many)`
    // let input_a_enc: NonInteractiveBatchedFheBools<Vec<Vec<u64>>> =

    // Server extracts the cyphertexts (1 cyphertext per bit)
    let now = std::time::Instant::now();
    let cts_a: Vec<_> = (0..DATE_ID_BITS + DATA_BITS)
        .map(|i| {
            let c = Encryptor::<_, NonInteractiveBatchedFheBools<Vec<Vec<u64>>>>::encrypt(
                ck_a,
                input_a.as_slice(),
            )
            .key_switch(0)
            .extract(i);
            FheBool { data: c }
        })
        .collect();
    let cts_b: Vec<_> = (0..DATE_ID_BITS + DATA_BITS)
        .map(|i| {
            let c = Encryptor::<_, NonInteractiveBatchedFheBools<Vec<Vec<u64>>>>::encrypt(
                ck_b,
                input_b.as_slice(),
            )
            .key_switch(1)
            .extract(i);
            FheBool { data: c }
        })
        .collect();
    println!("Server cyphertext extract time: {:?}", now.elapsed());

    // After extracting each client's private inputs, server proceeds to evaluate
    // the function
    let (cts_date_id_a, cts_data_a) = cts_a.as_slice().split_at(DATE_ID_BITS);
    let (cts_date_id_b, cts_data_b) = cts_b.as_slice().split_at(DATE_ID_BITS);
    let now = std::time::Instant::now();
    let cts_out = birthday_match_fhe_par(
        cts_date_id_a.try_into().unwrap(),
        cts_data_a.try_into().unwrap(),
        cts_date_id_b.try_into().unwrap(),
        cts_data_b.try_into().unwrap(),
    );
    println!("Birthday match FHE evaluation time: {:?}", now.elapsed());

    // Server has finished running compute. Clients can proceed to decrypt the
    // output ciphertext using multi-party decryption.

    // Client side //

    // In multi-party decryption, each client needs to come online, download output
    // ciphertext from the server, produce "output ciphertext" dependent decryption
    // share, and send it to other parties (either via p2p or via server). After
    // receving decryption shares from other parties, clients can independently
    // decrypt output ciphertext.

    // each client produces decryption share
    let now = std::time::Instant::now();
    let decryption_shares_0 = cts_out.0.clone().map(|ct| {
        cks.iter()
            .map(|k| k.gen_decryption_share(&ct))
            .collect::<Vec<_>>()
    });
    let decryption_shares_1 = cts_out.1.clone().map(|ct| {
        cks.iter()
            .map(|k| k.gen_decryption_share(&ct))
            .collect::<Vec<_>>()
    });
    println!("Clients decryption share gen time: {:?}", now.elapsed());

    // With all decryption shares, clients can aggregate the shares and decrypt the
    // ciphertext
    let now = std::time::Instant::now();
    let out_0: Vec<_> = cts_out
        .0
        .iter()
        .zip(decryption_shares_0.iter())
        .map(|(ct, dec_shares)| ck_a.aggregate_decryption_shares(ct, dec_shares))
        .collect();
    let out_1: Vec<_> = cts_out
        .1
        .iter()
        .zip(decryption_shares_1.iter())
        .map(|(ct, dec_shares)| ck_a.aggregate_decryption_shares(ct, dec_shares))
        .collect();
    println!("Client decrypt time: {:?}", now.elapsed());

    println!("out_0: {:x?}", binary_to_u64(&out_0));
    println!("out_1: {:x?}", binary_to_u64(&out_1));
}
