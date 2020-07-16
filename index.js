/*
LBM depends on:
1. advection (particles move)
2. collision (particles hit each other and try to relax themselves)

df/dt + e . divf/x = - 1/tau (f-f_eq)

LBM algo:
1. select site (node)
2. advection of particles from neighborhoods (movement)
3. collision of new set of particles that has been advected
4. post-collision rearrangement
5. next site (node) selection
*/

const { performance } = require("perf_hooks");
const { sleep } = require("./helpers");

const CELL_TYPE = "D2Q9";

// Speed of sound (m/s)
const vs = 334;
// Spatial step in cm
const spatialStep = 1;
// omega should be < 2 for positive viscosity
let omega = 1.0;
let rho0 = 1.0;
const accelX = 0.015;
const accelY = 0;

const inletDensity = 1.05;
const outletDensity = 1;
const inletSpeed = 0.5;
const outletSpeed = 0.5;
// Relaxation
const tau = 0.6;

const MAX_STEPS = 3;

const grid = {
  x: 2,
  y: 2,
};

let w0 = 4. / 9.;
let w1 = 1. / 9.;
let w2 = 1. / 36.;

let density = new Float32Array(9);

// 0, N, NE, E, SE, S, SW, W, NW (clockwise)
let weights = new Float32Array([w0, w1, w2, w1, w2, w1, w2, w1, w2]);
let ex = new Float32Array([0, 0, 1, 1, 1, 0, -1, -1, -1]);
let ey = new Float32Array([0, 1, 1, 0, -1, -1, -1, 0, 1]);

function cellSettings(type) {
  switch (type) {
    case "D2Q9":
      return { dimension: 2, speeds: 9 };
    case "D3Q19":
      return { dimension: 3, speeds: 19 };
    case "D3Q27":
      return { dimension: 3, speeds: 27 };
    default:
      return { dimension: 2, speeds: 9 };
  }
}

function getNextPos(direction, ex, ey, x, y) {
  return x + ex[direction], y - ey[direction];
}

function reverseDir(dir) {
  switch (dir) {
    case dir.N:
      return dir.S;
    case dir.S:
      return dir.N;
    case dir.W:
      return dir.E;
    case dir.E:
      return dir.W;
    case dir.NE:
      return dir.SW;
    case dir.SE:
      return dir.NW;
    case dir.NW:
      return dir.SE;
    case dir.SW:
      return dir.NE;
    default:
      return dir.none;
  }
}

function reflectDir(dir) {
  switch (dir) {
    case dir.N:
      return dir.S;
    case dir.S:
      return dir.N;
    case dir.W:
      return dir.E;
    case dir.E:
      return dir.W;
    case dir.NE:
      return dir.SE;
    case dir.SE:
      return dir.NE;
    case dir.NW:
      return dir.SW;
    case dir.SW:
      return dir.NW;
  }
}

function getTotalDensity() {
  let tDensity = 0;

  for (let i = 0; i < 9; i++) {
    tDensity += density[i];
  }

  return tDensity;
}

function getVelocity() {
  let tDensity = 0;
  let vx = 0;
  let vy = 0;

  for (let i = 0; i < 9; i++) {
    tDensity += density[i];
    vx += ex[i] * density[i];
    vy += ey[i] * density[i];
  }

  if (tDensity < 1e-14) return [0, 0];

  return [vx / tDensity, vy / tDensity];
}

// Important parts - collision and streaming (getting to equilibrium)

// this can be optimized, I won't do that to have the code easy to understand
// accelX, accelY are here to let you add a 'force' (as for example gravity, or some force to move the fluid at an inlet)
function equilibrium(accelXTau, accelYTau) {
  let totalDensity = density[0];
  // momentum (u)
  let ux = ex[0] * density[0];
  let uy = ey[0] * density[0];
  //console.log(ex[0], density[0], ex[0] * density[0], ux)
  // eqs 7 and 8 from https://iainhaslam.com/monplace/lbm-theory/
  for (let i = 1; i < 9; i++) {
    totalDensity += density[i];
    ux += ex[i] * density[i];
    uy += ey[i] * density[i];
  }

  if (totalDensity > 0) {
    ux /= totalDensity;
    uy /= totalDensity;
  }

  ux += accelXTau;

  //uy += accelYTau; // we don't want to move in y direction for now

  let u2 = ux * ux + uy * uy;
  // c ("coefficients") values are chosen so that particle densities "f_i"
  // propagate to nearby lattice nodes in precisely one timestep
  let c1 = 3.0;
  let c2 = 4.5; // 9. / 2.
  let c3 = -1.5; // - 3. / 2.

  let feq = new Float32Array(9);

  for (let i = 0; i < 9; i++) {
    // e_k * u
    let term = ex[i] * ux + ey[i] * uy;
    feq[i] = weights[i] * totalDensity * (1 + c1 * term + c2 * term * term + c3 * u2);
  }

  return feq;
}

async function collision(eqDist, tau) {
  for (let i = 0; i < 9; i++) {
    density[i] -= (density[i] - eqDist[i]) / tau;
  }

  return density;
}

async function fastCollision(eqDist, omega) {
  for (let i = 0; i < 9; i++) {
    density[i] -= omega * (density[i] - eqDist[i]);
  }

  return density;
}

async function collideAndStream() {
  // TODO: dont hardcode
  const latticeRows = grid.x;
  const latticeCols = grid.y;
  const lrMinus1 = latticeRows - 1;
  const lcMinus1 = latticeCols - 1;

  let accelXTau = accelX * tau;
  let accelYTau = accelY * tau;

  for (let y = 0; y < latticeRows; y++) {
    for (let x = 0; x < latticeCols; x++) {
      // get distribution function into equilibrium
      let eqDist = await equilibrium(accelXTau, accelYTau);
      // stream and collide
      await collision(eqDist, omega);
    }
  }
}

async function simulate() {
  const tStart = performance.now();

  const Cell = cellSettings(CELL_TYPE);

  await sleep(1000);

  // init
  for (let i = 0; i < Cell.speeds; i++) {
    density[i] = weights[i];
  }

  let processed = 0;

  for (let step = 0; step < MAX_STEPS; step++) {
    await collideAndStream();
    processed += grid.x * grid.y;
  }

  const tEnd = performance.now() - tStart;
  console.log(`Processed ${processed} updates in ${MAX_STEPS} time steps.`);
  console.log(`${tEnd}ms`);

  console.log(`${processed} updates per seconds.`);
}

simulate();
