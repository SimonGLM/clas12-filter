# Migration Todo List
| Function Name | Status | Comments |
| :--- | :--- | :--- |
| `ele_charge_cut` | ✅ | Covered by generic `charge_cut` |
| `prot_charge_cut` | ✅ | Covered by generic `charge_cut` |
| `neutr_charge_cut` | ✅ | Covered by generic `charge_cut` |
| `pip_charge_cut` | ✅ | Covered by generic `charge_cut` |
| `pim_charge_cut` | ✅ | Covered by generic `charge_cut` |
| `Kp_charge_cut` | ✅ | Covered by generic `charge_cut` |
| `Km_charge_cut` | ✅ | Covered by generic `charge_cut` |
| `phot_charge_cut` | ✅ | Covered by generic `charge_cut` |
| `FT_eid_charge_cut` | ✅ | Covered by generic `charge_cut` |
| `FT_photid_charge_cut` | ✅ | Covered by generic `charge_cut` |
| `CD_prot_charge_cut` | ❔ | Declared, not Defined |
| `CD_neutr_charge_cut` | ❔ | Declared, not Defined |
| `CD_pip_charge_cut` | ❔ | Declared, not Defined |
| `CD_pim_charge_cut` | ❔ | Declared, not Defined |
| `CD_Kp_charge_cut` | ❔ | Declared, not Defined |
| `CD_Km_charge_cut` | ❔ | Declared, not Defined |
| `ele_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `prot_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `neutr_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `pip_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `pim_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `Kp_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `Km_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `phot_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `FT_eid_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `FT_photid_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `CD_prot_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `CD_neutr_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `CD_pip_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `CD_pim_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `CD_Kp_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `CD_Km_default_PID_cut` | ✅ | Covered by generic `PID_cut` |
| `basic_FTOF_cut` | ✅ | Equal to "is somewhere in FTOF?", I guess..? |
| `Track_Quality_cut` | ❔ | Declared, not Defined |
| `CC_nphe_cut` | ✅ | |
| `EC_outer_vs_EC_inner_cut` | ✅ | Only checks energy from ECAL |
| `Km_EC_outer_vs_EC_inner_cut` | ✅ | Same as `EC_outer_vs_EC_inner_cut` with `loose` tightness |
| `pim_EC_outer_vs_EC_inner_cut` | ✅ | Same as `EC_outer_vs_EC_inner_cut` with `loose` tightness |
| `phot_EC_outer_vs_EC_inner_cut` | ✅ | Similar to `EC_outer_vs_EC_inner_cut` (checks ECin+ECout, instead of PCAL)|
| `EC_sampling_fraction_cut` | ✅ | Some TODOs in code: (bools & RCDB)
| `EC_hit_position_fiducial_cut` | ❔ | Declared, not Defined |
| `EC_hit_position_fiducial_cut_homogeneous` | ✅ | |
| `phot_EC_hit_position_fiducial_cut` | ✅ | Identical to `EC_hit_position_fiducial_cut_homogeneous`, except for the 13.5 value |
| `DC_fiducial_cut_XY` | ❔ | Declared, not Defined |
| `DC_fiducial_cut_theta_phi` | ❔ | Declared, not Defined |
| `DC_fiducial_cut_edge` | ✅ | |
| `DC_z_vertex_cut` | 👨‍💻 WIP | Exists in `cuts.h` (not implemented), already implemented in clas12::zVertexFilter? |
| `prot_delta_vz_cut` | ✅ | Covered by generic `delta_vz_cut` |
| `neutr_delta_vz_cut` | ✅ | Covered by generic `delta_vz_cut` |
| `pip_delta_vz_cut` | ✅ | Covered by generic `delta_vz_cut` |
| `pim_delta_vz_cut` | ✅ | Covered by generic `delta_vz_cut` |
| `Kp_delta_vz_cut` | ✅ | Covered by generic `delta_vz_cut` |
| `Km_delta_vz_cut` | ✅ | Covered by generic `delta_vz_cut` |
| `CD_prot_delta_vz_cut` | ✅ | Covered by generic `delta_vz_cut`, but different cuts (dep. on PID) |
| `CD_neutr_delta_vz_cut` | ✅ | Covered by generic `delta_vz_cut`, but different cuts (dep. on PID) |
| `CD_pip_delta_vz_cut` | ✅ | Covered by generic `delta_vz_cut`, but different cuts (dep. on PID) |
| `CD_pim_delta_vz_cut` | ✅ | Covered by generic `delta_vz_cut`, but different cuts (dep. on PID) |
| `CD_Kp_delta_vz_cut` | ✅ | Covered by generic `delta_vz_cut`, but different cuts (dep. on PID) |
| `CD_Km_delta_vz_cut` | ✅ | Covered by generic `delta_vz_cut`, but different cuts (dep. on PID) |
| `neutr_beta_cut` | ✅ | very similar to phot_beta_cut |
| `neutr_delta_beta_cut` | ❔ | Declared, not Defined |
| `neutr_tofmass_cut` | ❔ | Declared, not Defined |
| `pim_ele_reject_cut` | 🟡 N/A | checks if any electron checks has failed, but... why? |
| `Km_ele_reject_cut` | 🟡 N/A | Same as `pim_ele_reject_cut` |
| `maximum_probability_cut` | ❔ | Declared, not Defined |
| `phot_beta_cut` | ✅ | |
| `phot_EC_sampling_fraction_cut` | ✅ | |
| `FT_eid_FTCAL_fiducial_cut` | ✅ | |
| `FT_eid_FTTRK_fiducial_cut` | ✅ | Same as `FT_eid_FTCAL_fiducial_cut` |
| `FT_eid_FTHODO_fiducial_cut` | ✅ | Same as `FT_eid_FTCAL_fiducial_cut` |
| `FT_eid_energy_vs_radius_cut` | ✅ | Always `true`... why? |
| `FT_photid_FTCAL_fiducial_cut` | ✅ | No idea what this logic checks |
| `FT_photid_beta_cut` | ✅ | Always `true`... why? |
| `CD_prot_beta_cut` | ❔ | Declared, not Defined |
| `CD_neutr_beta_cut` | ✅ | Same as `neutr_beta_cut` |
| `CD_pip_beta_cut` | ❔ | Declared, not Defined |
| `CD_pim_beta_cut` | ❔ | Declared, not Defined |
| `CD_Kp_beta_cut` | ❔ | Declared, not Defined |
| `CD_Km_beta_cut` | ❔ | Declared, not Defined |