extern crate threshold_secret_sharing as tss;
extern crate time;
extern crate primal;

macro_rules! timer {
    ($what:expr, $card:expr, $code:block) => ({

        let start = time::get_time();
        $code
            let end = time::get_time();
        let sec = (end.sec - start.sec) as f64 +
            (end.nsec - start.nsec) as f64 / 1_000_000_000.0;

        println!("Ran {} {} in {}ms, {}ms each.", $card, $what, sec*1000.0, sec/($card as f64) * 1000.0);
    })
}

struct Server<'a> {
    params: &'a Params,
    superusers_pubkeys: Vec<u64>,
    superusers_inboxes: Vec<Vec<Vec<u64>>>,
    output_shares:Vec<Vec<(u64, u64)>>,
    result: Option<Vec<u64>>,
}

impl<'a> Server<'a> {
    fn new(params:&Params) -> Server {
        Server {
            params: params,
            superusers_pubkeys: std::iter::repeat(0).take(params.number_of_superusers).collect(),
            superusers_inboxes: std::iter::repeat(vec!()).take(params.number_of_superusers).collect(),
            output_shares: vec!(),
            result: None
        }
    }

    fn register_superuser_pubkey(&mut self, identity: usize, key: u64) {
        self.superusers_pubkeys.insert(identity, key);
    }

    fn get_superusers_pubkeys(&self) -> &[u64] {
        &self.superusers_pubkeys
    }

    fn upload_input_shares(&mut self, su:usize, data:Vec<u64>) {
        self.superusers_inboxes[su].push(data);
    }

    fn download_input_shares(&self, su:usize) -> &[Vec<u64>] {
        &*self.superusers_inboxes[su]
    }

    fn upload_output_shares(&mut self, data:Vec<(u64,u64)>) {
        self.output_shares.push(data);
    }

    fn combine(&mut self) {
        let data = (0..self.params.number_of_inputs).map(|i| {
            let v:Vec<(u64,u64)> = self.output_shares.iter().take(self.params.threshold).map(|data| data[i]).collect();
            self.params.tss.decrypt(&*v)
        }).collect();
        self.result = Some(data);
    }
}

#[derive(Debug)]
struct User<'a> {
    params:&'a Params,
    identity: usize,
    secrets: Vec<u64>,
    shares: Vec<Vec<(u64,u64)>> // shares[secret_id] = (x,y)
}

impl<'a> User<'a> {
    fn new(params:&Params, identity:usize) -> User {
        User { params:&params, identity:identity, shares:vec!(), secrets: params.make_secrets() }
    }

    fn generate_shares(&mut self) {
        for i in self.secrets.iter() {
            self.shares.push(self.params.secret_to_points(*i));
        }
    }

    fn upload_input_shares(&self, server:&mut Server) {
        let keys:Vec<u64> = server.get_superusers_pubkeys().iter().cloned().collect();
        for su in 0..self.params.number_of_superusers {
            let data:Vec<u64> = self.shares.iter().map( |secret_share| Params::crypt(keys[su], secret_share[su].1) ).collect();
            server.upload_input_shares(su, data);
        }
    }
}

#[derive(Debug)]
struct SuperUser<'a> {
    user: User<'a>,
    pubkey: u64,
    privkey: u64,
    input_shares: Option<Vec<Vec<u64>>>,
    output_shares: Option<Vec<(u64, u64)>>
}

impl<'a> SuperUser<'a> {
    fn new(params: &Params, identity: usize) -> SuperUser {
        SuperUser {
            user: User::new(&params, identity),
            pubkey: identity as u64 + 1,
            privkey: identity as u64 + 1,
            input_shares: None,
            output_shares: None,
        }
    }

    fn register(&self, server: &mut Server) {
        server.register_superuser_pubkey(self.user.identity, self.pubkey);
    }

    fn download_input_shares(&mut self, server: &Server) {
        let shares = server.download_input_shares(self.user.identity);
        self.input_shares = Some(
            shares.iter().map(|data|
                              data.iter().map(|&d|
                                              Params::decrypt(self.pubkey,d)).collect::<Vec<u64>>()
                             ).collect())
    }

    fn generate_output_shares(&mut self) {
        let prime = self.user.params.prime;
        let output_shares:Vec<(u64,u64)> =
            (0..self.user.params.number_of_inputs).map(|i| {
                ((self.user.identity+1) as u64,
                self.input_shares.as_ref().unwrap().iter().map(|part| part[i]).fold(0,|a,b|(a+b)%prime))
            }).collect();
        self.output_shares = Some(output_shares);
    }

    fn upload_output_shares(&self, server:&mut Server) {
        server.upload_output_shares(self.output_shares.as_ref().unwrap().clone());
    }
}

#[derive(Debug)]
struct Params {
    number_of_inputs:usize,
    number_of_users:usize,
    number_of_superusers:usize,
    input_bound:usize,
    output_bound:usize,
    threshold:usize,
    prime:u64,
    tss:tss::ThresholdSecretSharing
}

impl Params {
    fn new( number_of_inputs:usize, number_of_users:usize,
            number_of_superusers:usize, input_bound:usize,
            threshold:usize) -> Params {

        assert!(number_of_superusers < number_of_users);
        let output_bound = input_bound*number_of_users;


        let prime:u64 = primal::Sieve::new(9*output_bound).primes_from(8*output_bound).next().unwrap() as u64;
        //let prime:u64 = primal::Primes::all().filter( |&p| p > 8*output_bound).next().unwrap() as u64;

        Params {
            number_of_inputs:number_of_inputs, number_of_users:number_of_users,
            number_of_superusers:number_of_superusers, input_bound:input_bound,
            output_bound:output_bound, threshold:threshold, prime:prime,
            tss: tss::ThresholdSecretSharing::new(number_of_superusers, threshold, prime)
        }
    }

    fn make_secrets(&self) -> Vec<u64> {
        (0..self.number_of_inputs).map(|i| i as u64).collect()
    }

    fn secret_to_points(&self, secret:u64) -> Vec<(u64,u64)> {
        self.tss.encrypt(secret)
    }

    fn crypt(pubkey:u64, value:u64) -> u64 {
        //value
        pubkey ^ value
    }

    fn decrypt(privkey:u64, value:u64) -> u64 {
        //value
        privkey ^ value
    }
}


fn main() {
    let params = Params::new(2, 1000, 500, 1000, 20);

    println!("picked prime: {}", params.prime);

    let mut server = Server::new(&params);

    let mut superusers: Vec<_> = (0..params.number_of_superusers)
        .map(|i| SuperUser::new(&params, i))
        .collect();

    let mut regular_users: Vec<_> =
        (params.number_of_superusers..params.number_of_users).map(|i| User::new(&params,i)).collect();

    timer!("superuser registration", superusers.len(), {
        for su in superusers.iter() {
            su.register(&mut server);
        }
    });

    timer!("user share generation", params.number_of_users, {
        for mut u in regular_users.iter_mut() {
            u.generate_shares();
        }
        for mut u in superusers.iter_mut() {
            u.user.generate_shares();
        }
    });

    timer!("upload shares", params.number_of_users, {
        for u in regular_users.iter() {
            u.upload_input_shares(&mut server);
        }
        for u in superusers.iter() {
            u.user.upload_input_shares(&mut server);
        }
    });

    timer!("download shares", params.number_of_superusers, {
        for u in superusers.iter_mut() {
            u.download_input_shares(&mut server);
        }
    });

    timer!("sum shares", params.number_of_superusers, {
        for u in superusers.iter_mut() {
            u.generate_output_shares();
        }
    });

    timer!("upload shares", params.number_of_users, {
        for u in superusers.iter() {
            u.upload_output_shares(&mut server);
        }
    });

    timer!("server combine", 1, { server.combine() });

    println!("server final: {:?}", server.result.unwrap());
}

