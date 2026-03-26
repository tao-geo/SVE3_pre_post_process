# using ffmpeg to create video

# fig_prefix=case_A_ICE0_V1D_Vr_slice_lon270
# fig_prefix=case_B_ICE0_Vr_slice_lon270
fig_prefix=case_L17_ICE6G_V1D_Vr_slice_lon270
ffmpeg -framerate 1 -start_number 0 -i ../figs/${fig_prefix}.%d.png -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p ../figs/video/${fig_prefix}.mp4

# ffmpeg -framerate 3 -start_number 0 -i ./ANU_1x1_preprocessed/ice_final_minus_orig_%d.png -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p ANU_change_after_preprocessed.mp4

