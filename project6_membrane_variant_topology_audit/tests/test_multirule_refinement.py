"""Test rule-feature construction and constrained calibration with synthetic inputs."""
from __future__ import annotations
from copy import deepcopy
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"tools"))
from membrane_rule_audit.refinement import (ALPHABET,RULES,FeatureSpec,prepare_record,
    variant_features,apply_refinement,ablation_coefficients,fit_coefficients,topology_map)


def feature(kind,start,end=None,description="",direct=True):
    return dict(type=kind,location=dict(start=dict(value=start,modifier="EXACT"),
        end=dict(value=start if end is None else end,modifier="EXACT")),description=description,
        evidences=[dict(evidenceCode="ECO:0000269" if direct else "ECO:0000250",source="PubMed",id="synthetic")])


def fixture():
    seq = list("A"*45)
    seq[3:6] = list("NAS")
    seq[7] = "C"
    seq[11] = "C"
    record = dict(sequence="".join(seq),tmds=[dict(start=15,end=35)],
        loops={"1":dict(start=1,end=14,location="outside"),"2":dict(start=36,end=45,location="inside")})
    annotation = dict(primaryAccession="TEST",sequence=dict(value=record["sequence"]),
        features=[feature("Glycosylation",4,description="N-linked asparagine"),feature("Disulfide bond",8,12)])
    return record,annotation


def prepared():
    return prepare_record("TEST",*fixture())


def values(record,variants):
    return variant_features(record,[p for p,n,a in variants],[n for p,n,a in variants],[a for p,n,a in variants])


class TestAnnotations(unittest.TestCase):
    def test_three_additions_have_eligible_sources_without_free_controls(self):
        rec = prepared()
        self.assertEqual(len(rec["glycosites"]),1)
        self.assertEqual(len(rec["disulfides"]),1)
        self.assertFalse(rec["free_cysteine_controls_required"])

    def test_sequence_and_accession_mismatch_rejected(self):
        top,ann = fixture()
        with self.assertRaises(ValueError): prepare_record("OTHER",top,ann)
        ann["sequence"]["value"] += "A"
        with self.assertRaises(ValueError): prepare_record("TEST",top,ann)

    def test_inexact_bond_not_accepted(self):
        top,ann = fixture()
        ann["features"][1]["location"]["start"]["modifier"] = "UNSURE"
        self.assertFalse(prepare_record("TEST",top,ann)["disulfides"])

    def test_direct_and_all_annotated_arms_separate(self):
        top,ann = fixture()
        for f in ann["features"]: f["evidences"][0]["evidenceCode"] = "ECO:0000250"
        self.assertFalse(prepare_record("TEST",top,ann)["disulfides"])
        rec = prepare_record("TEST",top,ann,spec=FeatureSpec(annotation_evidence="all_annotated"))
        self.assertEqual(len(rec["disulfides"]),1)
        self.assertFalse(rec["disulfides"][0]["direct_evidence"])

    def test_cytosolic_motif_does_not_receive_glyco_feature(self):
        top,ann = fixture()
        top["loops"]["1"]["location"] = "inside"
        rec = prepare_record("TEST",top,ann)
        self.assertFalse(rec["glycosites"])
        self.assertEqual(values(rec,[(4,"N","K")])[0,2],0)

    def test_potential_unannotated_sequon_is_not_occupied(self):
        top,ann = fixture()
        ann["features"] = ann["features"][1:]
        self.assertEqual(values(prepare_record("TEST",top,ann),[(4,"N","A")])[0,2],0)

    def test_conflicting_partners_and_signal_sites_excluded(self):
        top,ann = fixture()
        seq = list(top["sequence"]); seq[9] = "C"; top["sequence"] = "".join(seq)
        ann["sequence"]["value"] = top["sequence"]
        ann["features"].append(feature("Disulfide bond",8,10))
        self.assertFalse(prepare_record("TEST",top,ann)["disulfides"])
        top,ann = fixture(); ann["features"].append(feature("Signal",1,12))
        rec = prepare_record("TEST",top,ann)
        self.assertFalse(rec["glycosites"]); self.assertFalse(rec["disulfides"])

    def test_binding_and_redox_excluded(self):
        top,ann = fixture(); ann["features"].append(feature("Binding site",8))
        self.assertFalse(prepare_record("TEST",top,ann)["disulfides"])
        top,ann = fixture(); ann["features"][1]["description"] = "redox switch"
        self.assertFalse(prepare_record("TEST",top,ann)["disulfides"])


class TestFeatures(unittest.TestCase):
    def test_topology_uses_training_encoding(self):
        geo = topology_map(fixture()[0])
        self.assertEqual([geo[p]["side"] for p in (1,15,40)],[2,0,1])

    def test_unknown_side_not_inside(self):
        top,ann = fixture(); top["loops"]["1"]["location"] = "unknown"
        rec = prepare_record("TEST",top,ann)
        self.assertEqual(values(rec,[(1,"A","K")])[0,0],0)

    def test_inside_outside_decay_and_scope(self):
        x = values(prepared(),[(14,"A","K"),(36,"A","K"),(14,"A","E"),(37,"A","R")])
        np.testing.assert_allclose(x[:,0],[.5,-.5,0,-.5*np.exp(-1/8)])

    def test_charge_depth_and_native_gates(self):
        x = values(prepared(),[(15,"A","K"),(20,"A","K"),(25,"A","K"),(25,"A","L"),(40,"A","K")])
        np.testing.assert_allclose(x[:,1],[0,.5,1,0,0])

    def test_native_charge_and_histidine_not_charge_gains(self):
        top,ann = fixture(); seq = list(top["sequence"]); seq[24] = "K"; seq[25] = "H"
        top["sequence"] = "".join(seq); ann["sequence"]["value"] = top["sequence"]
        x = values(prepare_record("TEST",top,ann),[(25,"K","R"),(26,"H","R"),(24,"A","H")])
        np.testing.assert_array_equal(x[:,1],0)

    def test_short_tm_not_forced(self):
        top,ann = fixture(); rec = prepare_record("TEST",top,ann,spec=FeatureSpec(minimum_tm_length=22))
        self.assertEqual(values(rec,[(25,"A","K")])[0,1],0)

    def test_glyco_three_positions_and_preserving_ST(self):
        x = values(prepared(),[(4,"N","A"),(5,"A","P"),(5,"A","G"),(6,"S","T"),(6,"S","A")])
        np.testing.assert_array_equal(x[:,2],[1,1,0,0,1])

    def test_all_disulfide_losses_but_not_free_cysteine_assumed(self):
        x = values(prepared(),[(8,"C","S"),(12,"C","A"),(8,"C","K"),(1,"A","C")])
        np.testing.assert_array_equal(x[:,3],[1,1,1,0])

    def test_native_identity_positions_and_hash_guarded(self):
        rec = prepared()
        for variants in ([(8,"A","S")],[(8,"C","C")],[(0,"A","C")],[(8.5,"C","A")]):
            with self.assertRaises(ValueError): values(rec,variants)
        rec["sequence_sha256"] = "bad"
        with self.assertRaises(ValueError): values(rec,[(8,"C","S")])


class TestCombination(unittest.TestCase):
    def test_additive_scope_and_zero_preservation(self):
        rec = prepared(); x = values(rec,[(4,"N","K"),(25,"A","R"),(8,"C","S"),(1,"A","L")])
        beta = dict(zip(RULES,[.2,.3,.4,.5])); raw = np.array([.2,.2,.2,.2])
        out = apply_refinement(raw,x,beta)
        np.testing.assert_allclose(out["scores"],raw+x@np.array(list(beta.values())))
        self.assertEqual(out["scores"][-1],raw[-1])
        np.testing.assert_array_equal(apply_refinement(raw,x,{r:0 for r in RULES})["scores"],raw)

    def test_no_clipping_nan_preserved_and_double_apply_rejected(self):
        x = np.ones((2,4)); beta = {r:1 for r in RULES}
        out = apply_refinement([.9,np.nan],x,beta)
        self.assertAlmostEqual(out["scores"][0],4.9); self.assertTrue(np.isnan(out["scores"][1]))
        with self.assertRaises(ValueError): apply_refinement([0,0],x,beta,already_refined=True)
        beta[RULES[0]] = -1
        with self.assertRaises(ValueError): apply_refinement([0,0],x,beta)

    def test_all_nine_declared_ablations(self):
        variants = ablation_coefficients({r:1 for r in RULES})
        self.assertEqual(len(variants),9)
        self.assertEqual(sum(variants["raw"].values()),0)
        self.assertTrue(all(v["positive_inside"] == 1 for k,v in variants.items() if k != "raw"))

    def test_disulfide_offset_does_not_fake_pairing_similarity(self):
        x = values(prepared(),[(8,"C","S"),(12,"C","S")]); raw = np.array([.1,.8])
        out = apply_refinement(raw,x,{r:.5 if r == "disulfide_loss" else 0 for r in RULES})
        self.assertAlmostEqual(abs(np.diff(out["scores"])[0]),abs(np.diff(raw)[0]))

    def test_all_alt_bond_loss_shift_preserves_within_position_contrast(self):
        alt = [a for a in ALPHABET if a != "C"]
        x = variant_features(prepared(),[8]*19,["C"]*19,alt)
        raw = np.linspace(.1,.9,19)
        post = apply_refinement(raw,x,{r:.3 if r == "disulfide_loss" else 0 for r in RULES})["scores"]
        kr = np.isin(alt,list("KR")); neutral = np.isin(alt,list("AFGILMNPQSTVWY"))
        self.assertAlmostEqual(post[kr].mean()-post[neutral].mean(),raw[kr].mean()-raw[neutral].mean())


class TestCalibration(unittest.TestCase):
    def data(self):
        rng = np.random.default_rng(27); x = rng.uniform(0,1,(120,4))
        acc = np.repeat(["a","b","c"],40)
        base = np.full(120,.2)
        kwargs = dict(accession=acc,assay=np.repeat("synthetic",120),family=acc.copy(),
                      role=np.repeat("calibration",120),target_scale="synthetic_ranking",
                      label_provenance="declared_synthetic_unit_test")
        return x,base,kwargs

    def test_known_positive_repair_without_free_cys_controls(self):
        x,base,kwargs = self.data(); beta = np.array([.2,.3,.4,.5])
        result = fit_coefficients(x,base,base+x@beta,ridge=1e-8,**kwargs)
        np.testing.assert_allclose(list(result["coefficients"].values()),beta,atol=1e-6)
        self.assertIn("control_free",result["disulfide_design"])

    def test_opposite_experimental_direction_cannot_force_positive_weight(self):
        x,base,kwargs = self.data()
        result = fit_coefficients(x,base,base-x[:,3],**kwargs)
        self.assertEqual(result["coefficients"]["disulfide_loss"],0)

    def test_passing_rule_can_be_locked_at_zero(self):
        x,base,kwargs = self.data()
        result = fit_coefficients(x,base,base+x.sum(axis=1),locked_zero_rules=["charge_burial"],**kwargs)
        self.assertEqual(result["coefficients"]["charge_burial"],0)

    def test_evaluation_rows_cannot_fit(self):
        x,base,kwargs = self.data(); kwargs["role"][0] = "evaluation"
        with self.assertRaises(ValueError): fit_coefficients(x,base,base,**kwargs)

    def test_unidentified_rule_reported_instead_of_claimed_pass(self):
        x,base,kwargs = self.data(); x[:,3] = 1
        out = fit_coefficients(x,base,base+x[:,2],**kwargs)
        self.assertIn("disulfide_loss",out["unsupported_rules"])
        self.assertEqual(out["coefficients"]["disulfide_loss"],0)

    def test_assay_offset_does_not_manufacture_bond_penalty(self):
        x,base,kwargs = self.data()
        out = fit_coefficients(x,base,base+np.repeat([5,-3,2],40),**kwargs)
        np.testing.assert_allclose(list(out["coefficients"].values()),0,atol=1e-12)

    def test_missing_and_family_inconsistency_rejected(self):
        x,base,kwargs = self.data(); target = base.copy(); target[0] = np.nan
        with self.assertRaises(ValueError): fit_coefficients(x,base,target,**kwargs)
        kwargs["family"][0] = "b"
        with self.assertRaises(ValueError): fit_coefficients(x,base,base,**kwargs)


if __name__ == "__main__":
    unittest.main()
