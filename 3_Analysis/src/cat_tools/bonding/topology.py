"""Compatibility imports; implementation lives in cat_shared.bonding.topology."""

from cat_shared.bonding.topology import (
    SiteImage as SiteImage,
    Contact as Contact,
    Polyhedron as Polyhedron,
    StructureData as StructureData,
    TopologyConfig as TopologyConfig,
    TopologyResult as TopologyResult,
    sha256_file as sha256_file,
    read_poscar as read_poscar,
    displacement_cart as displacement_cart,
    contact_distance as contact_distance,
    canonical_contact_key as canonical_contact_key,
    build_polyhedra as build_polyhedra,
    enumerate_pair_contacts as enumerate_pair_contacts,
    enumerate_inter_xx as enumerate_inter_xx,
    classify_contact as classify_contact,
    format_cohp_between as format_cohp_between,
    load_manifest as load_manifest,
    write_manifest as write_manifest,
    topology_config_from_manifest as topology_config_from_manifest,
    validate_manifest_structure as validate_manifest_structure,
)

__all__ = ['SiteImage', 'Contact', 'Polyhedron', 'StructureData', 'TopologyConfig', 'TopologyResult', 'sha256_file', 'read_poscar', 'displacement_cart', 'contact_distance', 'canonical_contact_key', 'build_polyhedra', 'enumerate_pair_contacts', 'enumerate_inter_xx', 'classify_contact', 'format_cohp_between', 'load_manifest', 'write_manifest', 'topology_config_from_manifest', 'validate_manifest_structure']
