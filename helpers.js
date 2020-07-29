export async function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

export function encodeVelocity(velocityX, velocityY) {
  var ux = Math.floor(velocityX);
  var uy = Math.floor(velocityY);
  return {
    r: (ux & 0xFF00) >> 8,
    g: (ux & 0x00FF),
    b: (uy & 0xFF00) >> 8,
    a: (uy & 0x00FF),
  }
}

export function fluidNodesLength(cells) {
  const fluidNodes = Array.from(cells.values()).filter(cell => cell.type === "fluid");
  return fluidNodes && fluidNodes.length || 0
}

export function getFluidNodes(cells) {
  return Array.from(cells.values()).filter(cell => cell.type === "fluid");
}

export function solidNodesLength(cells) {
  const solidNodes = Array.from(cells.values()).filter(cell => cell.type === "solid");
  return solidNodes && solidNodes.length || 0
}

export function getSolidNodes(cells) {
  return Array.from(cells.values()).filter(cell => cell.type === "solid");
}

// module.exports = {
//   sleep
// }
