use std::f64::INFINITY;
use std::ops::{Add, Div, Index, IndexMut, Mul, Sub};

#[derive(Debug, Copy, Clone, Serialize)]
pub struct Vector3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Index<usize> for Vector3D {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Index out of bounds. Vector3D only accepts indexes from 0 to 2."),
        }
    }
}

impl IndexMut<usize> for Vector3D {
    fn index_mut(&mut self, index: usize) -> &mut f64 {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Index out of bounds. Vector3D only accepts indexes from 0 to 2."),
        }
    }
}

impl Add for Vector3D {
    type Output = Self;

    fn add(self, other: Vector3D) -> Self::Output {
        Vector3D {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for Vector3D {
    type Output = Self;

    fn sub(self, other: Vector3D) -> Self::Output {
        Vector3D {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Mul<f64> for Vector3D {
    type Output = Self;

    fn mul(self, other: f64) -> Self::Output {
        Vector3D {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl Div<f64> for Vector3D {
    type Output = Self;

    fn div(self, other: f64) -> Self::Output {
        Vector3D {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

impl Vector3D {
    fn new<A: Into<f64>>(x: A, y: A, z: A) -> Vector3D {
        Vector3D {
            x: x.into(),
            y: y.into(),
            z: z.into(),
        }
    }

    pub fn magnitude(&self) -> f64 {
        (self.x.powf(2.0) + self.y.powf(2.0) + self.z.powf(2.0)).sqrt()
    }

    fn dummy_multiply(&self, vector: Vector3D) -> Vector3D {
        Vector3D {
            x: self.x * vector.x,
            y: self.y * vector.y,
            z: self.z * vector.z,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Particle {
    pub position: Vector3D,
    pub velocity: Vector3D,
    pub mass: f64,
    pub last_contact: f64,
}

impl Particle {
    pub fn new(position: Vector3D, velocity: Vector3D, mass: f64) -> Particle {
        Particle {
            position,
            velocity,
            mass,
            last_contact: 0.0,
        }
    }

    fn update(&mut self, time_delta: f64) {
        self.position = self.position + self.velocity * time_delta;
        self.last_contact += time_delta;
    }
}

struct Box {
    half_side_length: f64,
    particles: Vec<Particle>,
}

impl Box {
    fn new(volume: f64) -> Box {
        Box {
            half_side_length: volume.cbrt() * 0.5,
            particles: Vec::new(),
        }
    }

    fn side_area(&self) -> f64 {
        (self.half_side_length * 2.0).powf(2.0)
    }

    // Returns the time of contact and the reflection array
    fn calc_contact(&self, particle: Particle) -> (f64, Vector3D) {
        let mut t = INFINITY;

        let mut mult = Vector3D { x: 1.0, y: 1.0, z: 1.0 };

        let mut reflection_ind = 0;

        for i in 0..=2 {
            if particle.velocity[i] == 0.0 {
                continue;
            } else {
                let new_t = (particle.velocity[i].signum() * self.half_side_length
                    - particle.position[i])
                    / particle.velocity[i];

                if new_t < t {
                    t = new_t;
                    reflection_ind = i;
                }
            }
        }

        mult[reflection_ind] *= -1.0;

        (t, mult)
    }

    fn add_particle(&mut self, particle: Particle) {
        self.particles.push(particle);
    }
}

use serde_derive::Serialize;

#[derive(Clone, Copy, Debug, Serialize)]
pub struct Collision {
    pub position: Vector3D,
    pub wall_index: u8,
    pub impulse: f64,
    pub time: f64,
}

fn approx_equal(a: f64, b: f64) -> bool {
    (a.abs() - b.abs()).abs() < 2e-9
}

impl Collision {
    fn from_particle(
        part: &Particle,
        half_side_length: f64,
    ) -> Collision {
        let mut wall_index: u8 = 7;
        
        let reflection_sides = Vector3D {
            x: if approx_equal(part.position.x.abs(), half_side_length) {
                wall_index = if part.position.x.is_sign_negative() { 0 } else { 1 };
                1.0
            } else {
                0.0
            },
            y: if approx_equal(part.position.y.abs(), half_side_length) {
                wall_index = if part.position.y.is_sign_negative() { 2 } else { 3 };
                1.0
            } else {
                0.0
            },
            z: if approx_equal(part.position.z.abs(), half_side_length) {
                wall_index = if part.position.z.is_sign_negative() { 4 } else { 5 };
                1.0
            } else {
                0.0
            },
        };

        if wall_index == 7 {
            println!("HSL: {}, {:?}", half_side_length, part.position);
            wall_index = 0;
        }
        
        let reflection_vector = part.velocity.dummy_multiply(reflection_sides);        

        Collision {
            position: part.position,
            wall_index,
            impulse: 2.0 * reflection_vector.magnitude() * part.mass,
            time: part.last_contact,
        }
    }
}

use rand::{thread_rng, Rng};
use serde::Serialize;

// Boltzmann's constant
const K_B: f64 = 1.38064852e-23;

pub fn generate_position(half_side_length: f64) -> Vector3D {
    let x = thread_rng().gen_range(-half_side_length..=half_side_length);
    let y = thread_rng().gen_range(-half_side_length..=half_side_length);
    let z = thread_rng().gen_range(-half_side_length..=half_side_length);

    Vector3D { x, y, z }
}

use rand::prelude::*;
use rand_distr::StandardNormal;

pub fn generate_velocity(temperature: f64, mass: f64) -> Vector3D {
    let x = ((temperature * K_B) / mass).sqrt() * thread_rng().sample::<f64, _>(StandardNormal);
    let y = ((temperature * K_B) / mass).sqrt() * thread_rng().sample::<f64, _>(StandardNormal);
    let z = ((temperature * K_B) / mass).sqrt() * thread_rng().sample::<f64, _>(StandardNormal);

    Vector3D { x, y, z }
}

pub fn simulate (temp: f64, volume: f64, molecule_mass: f64, molecules_count: usize) -> (Vec<Vec<Collision>>, Vec<f64>) {
    let b = Box::new(volume);

    let mut all_sides: Vec<Vec<Collision>> = Vec::new();

    for _ in 0..6 {
        all_sides.push(Vec::new());
    }

    for _ in 0..molecules_count {  
        
        let mut part = Particle::new(
            generate_position(b.half_side_length),
            generate_velocity(temp, molecule_mass),
            molecule_mass
        );

        while part.last_contact < 1.0 {
            let (time_delta, reflection) = b.calc_contact(part);
            part.update(time_delta);

            let col = Collision::from_particle(&part, b.half_side_length);
            all_sides[col.wall_index as usize].push(col);

            part.velocity = part.velocity.dummy_multiply(reflection);
        }
    }

    for collisions in &mut all_sides {
        collisions.sort_by(|k, j| k.time.partial_cmp(&j.time).unwrap());
    }

    let mut pressures = Vec::new();

    for collisions in &mut all_sides {
        let mut cumulative_impulse = 0f64;
        let mut total_time = collisions.last().unwrap().time;

        for col in &*collisions {
            cumulative_impulse += col.impulse;
        }
        
        pressures.push((cumulative_impulse/total_time)/b.side_area());
    }

    return (all_sides, pressures);   
}
