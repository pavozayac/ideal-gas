mod simulation;
use simulation::*;

fn main() {
    let molecules_count = 10000;
    let molecule_mass = 5.3e-26;
    let volume = 3.0;
    let temperature = 300.0;

    simulate(temperature, volume, molecule_mass, molecules_count, rms_speed);

}