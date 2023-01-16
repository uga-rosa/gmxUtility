#!/usr/bin/env -S deno run -A

import { parse } from "https://deno.land/std@0.168.0/flags/mod.ts";
import { readLines } from "https://deno.land/std@0.168.0/io/mod.ts";
import { join } from "https://deno.land/std@0.168.0/path/mod.ts";

const REGION = 0.700;

type ParsedArgs = {
  structure_file: string;
  peg_z_coord_file: string;
  peg_rg_file: string;
  is_cyclic: boolean;
};

const parseArgs = (): ParsedArgs => {
  const parsedArgs = parse(Deno.args);
  return {
    structure_file: parsedArgs.f || "md-run.gro",
    peg_z_coord_file: parsedArgs.z || "peg_z.xvg",
    peg_rg_file: parsedArgs.r || "rg_peg.xvg",
    is_cyclic: parsedArgs.c || false,
  };
};

const readStructureFile = async (
  fileGro: string,
): Promise<[ausZcoord: number[], boxZ: number]> => {
  console.log(`Reading ${fileGro}...`);
  const lines = (await Deno.readTextFile(fileGro)).trim().split("\n");
  const boxSize = lines.pop()
    ?.trim()
    .split(/\s+/)
    .flatMap((e) => parseFloat(e) || []);
  if (!boxSize || boxSize.length < 3) {
    return Promise.reject(
      "Invalid file format: No box size in last line or less than 3 elements.",
    );
  }
  const regexAUS = /AUS/;
  const ausZcoordString = lines
    .filter((line) => regexAUS.test(line))
    // Cut only z coordinate.
    // https://manual.gromacs.org/documentation/current/reference-manual/file-formats.html#gro
    .map((line) => line.slice(36, 44).trim());
  const ausZcoord = [...new Set(ausZcoordString)]
    .flatMap((z) => parseFloat(z) || [])
    // ascending-order
    .sort();
  console.log("Done");
  return [ausZcoord, boxSize[2]];
};

const calcAdsorped = async (
  filePegZ: string,
  ausZcoord: number[],
  boxZ: number,
): Promise<number[][]> => {
  console.log("Calculating adsorption of PEG on GOLD...");
  const adsorpedPegIDs: number[][] = [];
  const file = await Deno.open(join(Deno.cwd(), filePegZ));
  for await (const line of readLines(file)) {
    // Skip comment line
    if (line.startsWith("@") || line.startsWith("#")) {
      continue;
    }
    const z_coord = line.trim()
      .split(/\s+/)
      .map((e) => parseFloat(e));
    // The first element is step number.
    z_coord.shift();
    const adsorpedPegID = z_coord.flatMap((z, idx) =>
      distance(z, ausZcoord, boxZ) < REGION ? idx : []
    );
    adsorpedPegIDs.push(adsorpedPegID);
  }
  console.log("Done");
  return adsorpedPegIDs;
};

const distance = (z: number, aus: number[], boxZ: number): number => {
  z -= aus[0];
  z = z > 0 ? z : z + boxZ;
  return Math.min(z - (aus[1] - aus[0]), boxZ - z);
};

const readRgFile = async (
  fileRg: string,
): Promise<number[]> => {
  console.log(`Reading ${fileRg}...`);
  const rgs: number[] = [];
  const file = await Deno.open(join(Deno.cwd(), fileRg));
  for await (const line of readLines(file)) {
    // Skip comment line
    if (line.startsWith("@") || line.startsWith("#")) {
      continue;
    }
    // step, Rg, Rgx, Rgy, Rgz
    const rg = line.trim()
      .split(/\s+/);
    rgs.push(parseFloat(rg[1]));
  }
  console.log("Done");
  return rgs;
};

const calcHistogram = (data: number[], bins: number): [number[], number[]] => {
  const max = getMaxOfArray(data);
  const min = getMinOfArray(data);
  const width = (max - min) / (bins - 1);
  const histogram: number[] = new Array(bins).fill(0);
  for (const item of data) {
    histogram[Math.floor((item - min) / width)]++;
  }
  const classValues = new Array(bins)
    .fill(min)
    .map((e, i) => e + (i + 0.5) * width);
  return [histogram, classValues];
};

const getMaxOfArray = (numArray: number[]): number =>
  numArray.reduce((a, b) => a > b ? a : b);

const getMinOfArray = (numArray: number[]): number =>
  numArray.reduce((a, b) => a < b ? a : b);

const writeHistogram = async (
  histogram: number[],
  classValues: number[],
  name: string,
) => {
  console.log("Writing the result of histogram...");
  const buffer = histogram
    .reduce(
      (prev, hist, idx) =>
        prev + [
          classValues[idx].toString().padEnd(18),
          hist,
        ].join(" ") + "\n",
      "",
    );
  await Deno.writeTextFile(`rgHist${name}.dat`, buffer);
  console.log("Done");
};

const writeResult = async (adsorpedPegIDs: number[][], rgs: number[]) => {
  console.log("Writing the result of adsorption...");
  const buffer = adsorpedPegIDs
    .reduce(
      (prev, ID, idx) =>
        prev +
        [
          idx.toString().padEnd(6, " "),
          rgs[idx].toString().padEnd(8, " "),
          ...ID.map((id) => id.toString().padStart(2, " ")),
        ].join(" ") + "\n",
      "",
    );
  await Deno.writeTextFile("adsorpedPegIDs.dat", buffer);
  console.log("Done");
};

const main = async () => {
  const parsedArgs = parseArgs();
  const [ausZcoord, boxZ] = await readStructureFile(parsedArgs.structure_file);
  const adsorpedPegIDs = await calcAdsorped(
    parsedArgs.peg_z_coord_file,
    ausZcoord,
    boxZ,
  );
  const steps = adsorpedPegIDs.filter((ID) => ID.length > 0).length;
  const allSteps = adsorpedPegIDs.length;
  console.log(
    `Adsorption probability: ${steps} / ${allSteps} = ${
      steps / allSteps * 100
    } %`,
  );

  const rgs = await readRgFile(parsedArgs.peg_rg_file);

  const [histogram, classValues] = calcHistogram(rgs, 50);
  await writeHistogram(histogram, classValues, "All");

  const adsorpedRgs = rgs
    .flatMap((e, i) => adsorpedPegIDs[i].length > 0 ? e : []);
  const [adsHistogram, adsClassValues] = calcHistogram(adsorpedRgs, 50);
  await writeHistogram(adsHistogram, adsClassValues, "Adsorped");

  const nonadsorpedRgs = rgs
    .flatMap((e, i) => adsorpedPegIDs[i].length === 0 ? e : []);
  const [nonadsHistogram, nonadsClassValues] = calcHistogram(
    nonadsorpedRgs,
    50,
  );
  await writeHistogram(nonadsHistogram, nonadsClassValues, "Nondsorped");

  await writeResult(adsorpedPegIDs, rgs);
};

main();
