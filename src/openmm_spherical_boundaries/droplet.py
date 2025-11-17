"""Utility for building a standalone spherical water droplet."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

from openmm import Vec3, XmlSerializer, unit
from openmm.app import CutoffNonPeriodic, Element, ForceField, Modeller, PDBFile, Topology

from .utils import clean_openmm_positions, distance

LOGGER = logging.getLogger(__name__)

DEFAULT_FORCEFIELD = ("amber99sb.xml", "spce.xml")


def prepare_water_droplet(
    radius: float,
    output_pdb: Path | str = "droplet.pdb",
    system_xml: Path | str = "droplet.xml",
    forcefield_files: Sequence[str] | None = None,
    solvent_model: str = "spce",
    padding: float = 0.3,
) -> int:
    """Create a spherical water droplet by solvating a dummy system and trimming."""

    if radius <= 0:
        raise ValueError("Radius must be positive.")

    if forcefield_files is None:
        forcefield_files = DEFAULT_FORCEFIELD
    forcefield_tuple = tuple(forcefield_files)
    if set(forcefield_tuple) != set(DEFAULT_FORCEFIELD):
        raise NotImplementedError(
            "Droplet generation currently supports only amber99sb.xml + spce.xml."
        )

    forcefield = ForceField(*forcefield_tuple)
    topology, positions = _single_water()
    modeller = Modeller(topology, positions)

    box_edge = 2.0 * (radius + padding)
    LOGGER.info(
        "Packing water using %s model into %.2f nm cube before trimming", solvent_model, box_edge
    )
    modeller.addSolvent(
        forcefield,
        model=solvent_model,
        boxSize=Vec3(box_edge, box_edge, box_edge),
    )

    trimmed_count = _trim_to_radius(modeller, radius)
    LOGGER.info("Trimmed droplet contains %d water molecules", trimmed_count)

    system = forcefield.createSystem(
        modeller.getTopology(),
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=1 * unit.nanometer,
        removeCMMotion=True,
    )

    output_pdb_path = Path(output_pdb)
    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    with output_pdb_path.open("w") as handle:
        PDBFile.writeFile(modeller.getTopology(), modeller.getPositions(), handle, keepIds=True)

    system_xml_path = Path(system_xml)
    system_xml_path.parent.mkdir(parents=True, exist_ok=True)
    with system_xml_path.open("w") as handle:
        handle.write(XmlSerializer.serialize(system))

    return trimmed_count


def _single_water() -> tuple[Topology, list[Vec3]]:
    """Return a topology and coordinates for a single SPC/E-like water molecule."""

    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("HOH", chain)
    element_o = Element.getBySymbol("O")
    element_h = Element.getBySymbol("H")
    atom_o = topology.addAtom("O", element_o, residue)
    atom_h1 = topology.addAtom("H1", element_h, residue)
    atom_h2 = topology.addAtom("H2", element_h, residue)
    topology.addBond(atom_o, atom_h1)
    topology.addBond(atom_o, atom_h2)

    positions = [
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.09572, 0.0, 0.0),
        Vec3(-0.031906, 0.092386, 0.0),
    ]
    return topology, unit.Quantity(positions, unit.nanometer)


def _trim_to_radius(modeller: Modeller, radius: float) -> int:
    """Remove residues whose oxygen atoms lie outside the requested radius."""

    pos_array = clean_openmm_positions(modeller.getPositions())
    residues_to_remove = []
    for residue in modeller.getTopology().residues():
        oxygen_index = None
        for atom in residue.atoms():
            if getattr(atom.element, "symbol", None) == "O":
                oxygen_index = atom.index
                break
        if oxygen_index is None:
            continue
        coord = pos_array[oxygen_index]
        if distance(coord, (0.0, 0.0, 0.0)) > radius:
            residues_to_remove.append(residue)

    if residues_to_remove:
        LOGGER.info("Removing %d residues outside %.2f nm sphere", len(residues_to_remove), radius)
        modeller.delete(residues_to_remove)

    return modeller.getTopology().getNumResidues()
