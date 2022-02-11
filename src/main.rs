mod simulation;
use chrono::{DateTime, Utc};
use serde_derive::Serialize;
use simulation::*;
use std::{fs, io::BufWriter, path::Path};

#[derive(Clone, Copy, Debug, Serialize)]
pub struct CollisionRecord {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub wall_index: u8,
    pub impulse: f64,
    pub time: f64,
}

fn main() {
    let molecules_count = 10000;
    let molecule_mass = 5.3e-26;
    let volume = 0.001;
    let temperature = 300.0;

    let args: Vec<String> = std::env::args().collect();

    if args[1] == "gay-lussac" {
        let dir_path = &format!(
            "Gay_Lussac_{}.csv",
            Utc::now().to_rfc3339().replace(":", "-").replace("-", "_").replace("+", "_").replace(".", "_")
        );

        println!("{}", dir_path);

        fs::create_dir(Path::new(dir_path)).expect("Could not create dir.");
        let pressures_file = fs::File::create(Path::new(&format!(
            "{}/pressures.csv",
            dir_path
        )))
        .expect("Cannot open file");
        let mut pressure_bw = BufWriter::new(pressures_file);
        let mut csv_p = csv::Writer::from_writer(pressure_bw);

        csv_p
            .write_record(&[
                "Side 0",
                "Side 1",
                "Side 2",
                "Side 3",
                "Side 4",
                "Side 5",
                "Temperature",
                "Volume",
            ])
        .expect("Could not write record.");

        csv_p.flush().expect("No flush");

        let mut all_pressures = Vec::new();

        for temp in args[2].parse::<i32>().unwrap()..=args[3].parse::<i32>().unwrap() {
            let t = temp;
            let (collisions, mut pressures) = simulate(
                temp as f64,
                volume,
                molecule_mass,
                molecules_count,
                rms_speed,
            );
            pressures.push(temp as f64);
            pressures.push(volume as f64);

            all_pressures.push(pressures);

            if temp == 1 || temp == 500 {
                for side in collisions {
                    let ind = side[0].wall_index;
                    let collisions_file = fs::File::create(Path::new(&format!(
                        "{}/t-{}-s-{}-collisions.csv",
                        dir_path,
                        t,
                        side[0].wall_index
                    )))
                    .expect("Cannot open file");
                    let mut collisions_bw = BufWriter::new(collisions_file);
                    let mut csv_c = csv::Writer::from_writer(collisions_bw);

                    for col in side {
                        let record = CollisionRecord {
                            x: col.position.x,
                            y: col.position.y,
                            z: col.position.z,
                            impulse: col.impulse,
                            time: col.time,
                            wall_index: col.wall_index
                        };

                        csv_c.serialize(record).expect("Could not serialize");
                    }

                    csv_c.flush().expect("Could not flush.");
                }
            }
        }

        

        for pressures in all_pressures {
            csv_p.serialize(&pressures).expect("Could not serialize.");
        }

        csv_p.flush().expect("Could not flush.");
    } else if args[1] == "boyle" {
        let dir_path = &format!(
            "Boyle_{}.csv",
            Utc::now().to_rfc3339().replace(":", "-").replace("-", "_").replace("+", "_").replace(".", "_")
        );

        println!("{}", dir_path);

        fs::create_dir(Path::new(dir_path)).expect("Could not create dir.");
        let pressures_file = fs::File::create(Path::new(&format!(
            "{}/pressures.csv",
            dir_path
        )))
        .expect("Cannot open file");
        let mut pressure_bw = BufWriter::new(pressures_file);
        let mut csv_p = csv::Writer::from_writer(pressure_bw);

        csv_p
            .write_record(&[
                "Side 0",
                "Side 1",
                "Side 2",
                "Side 3",
                "Side 4",
                "Side 5",
                "Temperature",
                "Volume",
            ])
        .expect("Could not write record.");

        csv_p.flush().expect("No flush");

        let mut all_pressures = Vec::new();

        for vol in args[2].parse::<i32>().unwrap()..=args[3].parse::<i32>().unwrap() {
            let v = 0.001 * vol as f64;
            let (collisions, mut pressures) = simulate(
                temperature,
                v,
                molecule_mass,
                molecules_count,
                rms_speed,
            );
            pressures.push(temperature as f64);
            pressures.push(v as f64);

            all_pressures.push(pressures);

            if vol == 1 || vol == 500 {
                for side in collisions {
                    let ind = side[0].wall_index;
                    let collisions_file = fs::File::create(Path::new(&format!(
                        "{}/v-{}-s-{}-collisions.csv",
                        dir_path,
                        v,
                        side[0].wall_index
                    )))
                    .expect("Cannot open file");
                    let mut collisions_bw = BufWriter::new(collisions_file);
                    let mut csv_c = csv::Writer::from_writer(collisions_bw);
    
                    for col in side {
                        let record = CollisionRecord {
                            x: col.position.x,
                            y: col.position.y,
                            z: col.position.z,
                            impulse: col.impulse,
                            time: col.time,
                            wall_index: col.wall_index
                        };
    
                        csv_c.serialize(record).expect("Could not serialize");
                    }
    
                    csv_c.flush().expect("Could not flush.");
                }            }
        }

        

        for pressures in all_pressures {
            csv_p.serialize(&pressures).expect("Could not serialize.");
        }

        csv_p.flush().expect("Could not flush.");
    } else {
        panic!("No experiment chosen.");
    }
}
