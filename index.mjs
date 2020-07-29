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

import {
  sleep,
  fluidNodesLength,
  getFluidNodes,
  solidNodesLength,
  getSolidNodes
} from "./helpers.js";

const CELL_TYPE = "D2Q9";

const LatticeCells = new Map();

// Speed of sound (m/s)
const vs = 334;
// Speed of sound squared
const cSq = Math.pow(vs, 2);
// Spatial step in cm
const spatialStep = 1;
// omega should be < 2 for positive viscosity
let omega = 1.0;
let rho0 = 1.0;
const accelX = 1e-6; // 0.015;
const accelY = 0;

const inletDensity = 1.05;
const outletDensity = 1;
const inletSpeed = 0.5;
const outletSpeed = 0.5;
// Relaxation
const tau = 0.6;

const T_MAX = 5;

const STEPS_IN_FRAME = 5;
let isGhost = false;

const boxSize = 20; // px
const rows = 20;
const cols = 60;

// TODO: dont use this for now
// const grid = {
//   x: 2,
//   y: 2,
// };

let w0 = 4 / 9; // 0 (center)
let w1 = 1 / 9; // N, E, S, W
let w2 = 1 / 36; // NE, SE, NW, SW (diagonal)

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

function initLattice() {
  let Lattice = {
    cells: new Map(),
    solidNodesIndices: [],
    nx: cols,
    ny: rows
  };

  for (var y = 0; y < Lattice.ny; y++) {
    for (var x = 0; x < Lattice.nx; x++) {
      Lattice.cells.set(`${x},${y}`, {
        x,
        y,
        type: "fluid",
        velocity: [0, 0],
        ex: new Float32Array([0, 0, 1, 1, 1, 0, -1, -1, -1]),
        ey: new Float32Array([0, 1, 1, 0, -1, -1, -1, 0, 1]),
        density: new Float32Array(weights) // initialize with initial weights added in
      });
    }
  }

  // solid
  for (let i = 0; i < 100; i++) {
    const x = Math.floor(Math.random() * cols);
    const y = Math.floor(Math.random() * rows);

    Lattice.cells.set(`${x},${y}`, {
      x,
      y,
      type: "solid",
      velocity: [0, 0]
    })

    Lattice.solidNodesIndices.push([x, y]);
  }

  return Lattice;
}

const LATTICE = initLattice();
const GHOST_LATTICE = Object.assign({}, LATTICE);

const CellDirection = {
  none: 0,
  N: 1,
  NE: 2,
  E: 3,
  SE: 4,
  S: 5,
  SW: 6,
  W: 7,
  NW: 8
}

const directions = Object.entries(CellDirection);

// function getNextPos(directionIndex, x, y) {
//   if (cell.type === "solid") return;
//   return [x + cell.ex[directionIndex], cell.y + cell.ey[directionIndex]];
// }

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

function getDensity(densities) {
  let tDensity = 0;

  for (let i = 0; i < 9; i++) {
    tDensity += densities[i];
  }

  return tDensity;
}

function getVelocity(densities) {
  let tDensity = 0;
  let vx = 0;
  let vy = 0;

  for (let i = 0; i < 9; i++) {
    tDensity += densities[i];
    vx += ex[i] * densities[i];
    vy += ey[i] * densities[i];
  }

  if (tDensity < 1e-14) return [0, 0];

  return [vx / tDensity, vy / tDensity];
}

// alternative way of getting the velocity
// https://github.com/mrjleo/LBM-js/blob/master/src/js/lbm/lbm.js#L155-L159
// function getVelocity(density, cell) {
//   // TODO: Vec2 impl. should go outside
//   function Vec2(x, y) {
//     this.x = x;
//     this.y = y;

//     this.mult = function(val) {
//       this.x *= val;
//       this.y *= val;
//     };
//   }

// 	let velocity = new Vec2(cell[2] + cell[3] + cell[4] - cell[8] - cell[7] - cell[6], cell[6] + cell[5] + cell[4] - cell[8] - cell[1] - cell[2]);
//   // total density in a given cell, should be calculated beforehand
//   velocity.mult(1 / density);
// 	return velocity; // { x: number, y: number }
// }

// this can be optimized, I won't do that to have the code easy to understand
// accelX, accelY are here to let you add a 'force' (as for example gravity, or some force to move the fluid at an inlet)
async function equilibrium(densities, accelXTau, accelYTau) {
  let totalDensity = densities[0];
  // momentum (u)
  let ux = ex[0] * densities[0];
  let uy = ey[0] * densities[0];
  //console.log(ex[0], density[0], ex[0] * density[0], ux)
  // eqs 7 and 8 from https://iainhaslam.com/monplace/lbm-theory/
  for (let i = 1; i < 9; i++) {
    totalDensity += densities[i];
    ux += ex[i] * densities[i];
    uy += ey[i] * densities[i];
  }

  if (totalDensity > 0) {
    ux /= totalDensity;
    uy /= totalDensity;
  }

  if (accelXTau) {
    ux += accelXTau;
  }

  // we don't want to move in y direction for now
  // if (accelXTau) {
  //   uy += accelYTau;
  // }

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

function getNextPos(directionIndex, x, y) {
  return `${x + ex[directionIndex]},${y - ex[directionIndex]}`
}

// main and ghost lattice
async function propagate(latticeSrc, latticeDest) {
  for (let x = 1; x < latticeSrc.length - 1; x++) {
    for (let y = 1; y < latticeSrc.length - 1; y++) {
      const oldSrcCell = latticeSrc.cells.get(`${x},${y}`);
      if (oldSrcCell.type !== 'solid') {
        const newDensities = [
          latticeSrc.cells.get(getNextPos(0)).density[0], // 0
          latticeSrc.cells.get(getNextPos(1)).density[1], // N
          latticeSrc.cells.get(getNextPos(2)).density[2], // NE
          latticeSrc.cells.get(getNextPos(3)).density[3], // E
          latticeSrc.cells.get(getNextPos(4)).density[4], // SE
          latticeSrc.cells.get(getNextPos(5)).density[5], // S
          latticeSrc.cells.get(getNextPos(6)).density[6], // SW
          latticeSrc.cells.get(getNextPos(7)).density[7], // W
          latticeSrc.cells.get(getNextPos(8)).density[8], // NW
        ];

        latticeDest.cells.set(`${x},${y}`, {
          ...oldSrcCells,
          density: newDensities
        })
      }
    }
  }

  bouncebackObstacles(latticeSrc, latticeDest);
}

function bouncebackObstacles(latticeSrc, latticeDest) {
  for (let x = 1; x < latticeSrc.length - 1; x++) {
    for (let y = 1; y < latticeSrc.length - 1; y++) {
      const oldSrcCell = latticeSrc.cells.get(`${x},${y}`);
      if (oldSrcCell.type === 'solid') {
        // Reverse directions
        const newDensities = [
          latticeSrc.cells.get(getNextPos(0)).density[0], // 0
          latticeSrc.cells.get(getNextPos(1)).density[5], // N -> S
          latticeSrc.cells.get(getNextPos(2)).density[6], // NE - SW
          latticeSrc.cells.get(getNextPos(3)).density[7], // E -> W
          latticeSrc.cells.get(getNextPos(4)).density[8], // SE - NW
          latticeSrc.cells.get(getNextPos(5)).density[1], // S -> N
          latticeSrc.cells.get(getNextPos(6)).density[2], // SW -> NE
          latticeSrc.cells.get(getNextPos(7)).density[3], // W -> E
          latticeSrc.cells.get(getNextPos(8)).density[4], // NW -> SE
        ];

        latticeDest.cells.set(`${x},${y}`, {
          ...oldSrcCells,
          density: newDensities
        })
      }
    }
  }
}

async function collision(densities, eqDist, tau) {
  for (let i = 0; i < 9; i++) {
    densities[i] -= (densities[i] - eqDist[i]) / tau;
  }

  return densities;
}

async function fastCollision(densities, eqDist, omega) {
  for (let i = 0; i < 9; i++) {
    densities[i] -= omega * (densities[i] - eqDist[i]);
  }

  return densities;
}

async function collideAndStream(lattice, ghostLattice) {
  const latticeRows = rows;
  const latticeCols = cols;

  let accelXTau = accelX * tau;
  let accelYTau = accelY * tau;

  let currentLattice = lattice;

  for (let y = 0; y < latticeRows; y++) {
    const shouldCollide = y !== 0 && y !== latticeRows - 1;

    for (let x = 0; x < latticeCols; x++) {
      for (var i = 0; i < STEPS_IN_FRAME; i++) {
        if (isGhost) {
          await propagate(lattice, ghostLattice);
          currentLattice = ghostLattice;
        } else {
          await propagate(ghostLattice, lattice);
          currentLattice = lattice;
        }
        isGhost = !isGhost;

        const cell = currentLattice.cells.get(`${x},${y}`);
        if (cell && cell.type === "fluid"/* && shouldCollide*/) { // TODO: use shouldCollide when other bouncebacks are implemented
          // get distribution function into equilibrium
          let eqDist = await equilibrium(cell.density);

          // collide
          const densitiesAfterCollision = await collision(cell.density, eqDist, omega);
          const velocityVec = getVelocity(densitiesAfterCollision); // [u, v]

          lattice.cells.set(`${x},${y}`, { ...cell, density: densitiesAfterCollision, velocity: velocityVec });
        }
      }

      // await sleep(500);
    }
  }

  return currentLattice
}

/* SIMULATION CODE */

async function* simulate() {
  const tStart = performance.now();

  await sleep(1000);

  // let averageVelocity = 1;
  // let prevAverageVelocity = 1;
  // const fluidNodes = getFluidNodes(lattice.cells).length;

  let processed = 0;

  for (let t = 0; t < T_MAX; t++) {
    await collideAndStream(LATTICE, GHOST_LATTICE);
    processed += rows * cols;

    // prevAverageVelocity = averageVelocity;
    // averageVelocity = getFluidNodes(lattice.cells).

    // yield LATTICE;
  }

  yield LATTICE;
  console.log(LATTICE)

  const tEnd = performance.now() - tStart;
  console.log(`Processed ${processed} updates in ${T_MAX} time steps.`);
  console.log(`${tEnd}ms`);

  console.log(`${processed} updates per seconds.`);
}

export { simulate, LATTICE }
