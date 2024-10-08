# Makefile generated by minimake

all: FICTURE/minimal_nthread8/batched.matrix.tsv.gz FICTURE/minimal_nthread8/coordinate_minmax.tsv FICTURE/minimal_nthread8/hexagon.d_15.tsv.gz FICTURE/minimal_nthread8/analysis/nF42.d_15/nF42.d_15.done FICTURE/minimal_nthread8/analysis/nF42.d_15/nF42.d_15.decode.prj_15.r_4_5.done

FICTURE/minimal_nthread8/batched.matrix.tsv.gz: FICTURE/MouseBrain.tsv.gz
	$(info --------------------------------------------------------------)
	$(info Creating minibatch from FICTURE/MouseBrain.tsv.gz...)
	$(info --------------------------------------------------------------)
	ficture make_spatial_minibatch --input FICTURE/MouseBrain.tsv.gz --output FICTURE/minimal_nthread8/batched.matrix.tsv --mu_scale 2.0 --batch_size 500 --batch_buff 30 --major_axis X
	sort -S 4G -k 2,2n -k 1,1g FICTURE/minimal_nthread8/batched.matrix.tsv | gzip -c > FICTURE/minimal_nthread8/batched.matrix.tsv.gz
	rm FICTURE/minimal_nthread8/batched.matrix.tsv

FICTURE/minimal_nthread8/coordinate_minmax.tsv: FICTURE/MouseBrain.tsv.gz
	$(info --------------------------------------------------------------)
	$(info Obtaining boundary coordinates to FICTURE/minimal_nthread8/coordinate_minmax.tsv...)
	$(info --------------------------------------------------------------)
	bash FICTURE/minimal_nthread8/write_minmax.sh FICTURE/MouseBrain.tsv.gz FICTURE/minimal_nthread8/coordinate_minmax.tsv 2.0

FICTURE/minimal_nthread8/hexagon.d_15.tsv.gz: FICTURE/MouseBrain.tsv.gz
	$(info --------------------------------------------------------------)
	$(info Creating DGE for 15um...)
	$(info --------------------------------------------------------------)
	ficture make_dge --key Count --input FICTURE/MouseBrain.tsv.gz --output FICTURE/minimal_nthread8/hexagon.d_15.tsv --hex_width 15 --n_move 2 --min_ct_per_unit 50 --mu_scale 2.0 --precision 2 --major_axis X
	sort -S 4G -k 1,1n FICTURE/minimal_nthread8/hexagon.d_15.tsv | gzip -c > FICTURE/minimal_nthread8/hexagon.d_15.tsv.gz
	rm FICTURE/minimal_nthread8/hexagon.d_15.tsv

FICTURE/minimal_nthread8/analysis/nF42.d_15/nF42.d_15.done: FICTURE/MouseBrain.tsv.gz FICTURE/minimal_nthread8/hexagon.d_15.tsv.gz
	$(info --------------------------------------------------------------)
	$(info Creating LDA for 15um and 42 factors...)
	$(info --------------------------------------------------------------)
	mkdir -p FICTURE/minimal_nthread8/analysis/nF42.d_15/figure
	ficture fit_model --input FICTURE/minimal_nthread8/hexagon.d_15.tsv.gz --output FICTURE/minimal_nthread8/analysis/nF42.d_15/nF42.d_15  --nFactor 42 --epoch 3 --epoch_id_length 2 --unit_attr X Y --key Count --min_ct_per_feature 20 --test_split 0.5 --R 10 --thread 8
	touch FICTURE/minimal_nthread8/analysis/nF42.d_15/nF42.d_15.done

FICTURE/minimal_nthread8/analysis/nF42.d_15/nF42.d_15.decode.prj_15.r_4_5.done: FICTURE/minimal_nthread8/batched.matrix.tsv.gz FICTURE/minimal_nthread8/hexagon.d_15.tsv.gz FICTURE/minimal_nthread8/analysis/nF42.d_15/nF42.d_15.done
	$(info --------------------------------------------------------------)
	$(info Creating projection for 15um and 42 factors, at 15um)
	$(info --------------------------------------------------------------)
	ficture transform --input FICTURE/MouseBrain.tsv.gz --output_pref FICTURE/minimal_nthread8/analysis/nF42.d_15/nF42.d_15.prj_15.r_4 --model FICTURE/minimal_nthread8/analysis/nF42.d_15/nF42.d_15.model.p --key Count --major_axis X --hex_width 15 --n_move 3 --min_ct_per_unit 20 --mu_scale 2.0 --thread 8 --precision 2
	$(info --------------------------------------------------------------)
	$(info Performing pixel-level decoding..)
	$(info --------------------------------------------------------------)
	ficture slda_decode --input FICTURE/minimal_nthread8/batched.matrix.tsv.gz --output FICTURE/minimal_nthread8/analysis/nF42.d_15/nF42.d_15.decode.prj_15.r_4_5 --model FICTURE/minimal_nthread8/analysis/nF42.d_15/nF42.d_15.model.p --anchor FICTURE/minimal_nthread8/analysis/nF42.d_15/nF42.d_15.prj_15.r_4.fit_result.tsv.gz --anchor_in_um --neighbor_radius 5 --mu_scale 2.0 --key Count --precision 0.1 --lite_topk_output_pixel 3 --lite_topk_output_anchor 3 --thread 8
	touch FICTURE/minimal_nthread8/analysis/nF42.d_15/nF42.d_15.decode.prj_15.r_4_5.done
