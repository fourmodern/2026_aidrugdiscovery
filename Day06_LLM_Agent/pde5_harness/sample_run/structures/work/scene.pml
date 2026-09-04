
load /home/hjpark/lecture_drug1/0905_agent/pde5_harness/sample_run/structures/work/receptor.pdb, rec
load /home/hjpark/lecture_drug1/0905_agent/pde5_harness/sample_run/structures/work/ligand_ref.sdf, refl
load /home/hjpark/lecture_drug1/0905_agent/pde5_harness/sample_run/structures/work_controlled/CHEMBL1916483_render_m1.sdf, dock

bg_color white
set ray_opaque_background, 1
set cartoon_transparency, 0.15
set stick_radius, 0.14
set sphere_scale, 0.28
set label_size, 16
set label_color, black
set antialias, 2
set ray_shadows, 0
set depth_cue, 0
set specular, 0.2

hide everything
select pocket, byres (rec within 5 of refl) and polymer
show sticks, pocket
color palecyan, pocket and elem C
show sticks, refl
color grey60, refl and elem C
show sticks, dock
color marine, dock and elem C
show spheres, rec and (resn ZN+MG)
color orange, rec and resn ZN
color green, rec and resn MG
distance polar, dock, pocket, 3.5, mode=2
color red, polar
set dash_width, 2.5
set dash_gap, 0.35
select shown, pocket and (resi 816+820+817+612+768+782+813+661+682+723+764+765)
label shown and name CB, "%s%s" % (resn, resi)
set label_size, 19
set label_outline_color, white
set label_position, (0, 0, 2.5)
hide sticks, pocket and not shown
orient refl
zoom refl, 3.2
turn y, 12

png sample_run/report/figures_controlled/fig8_binding_pocket.png, width=1800, height=1200, dpi=300, ray=1
