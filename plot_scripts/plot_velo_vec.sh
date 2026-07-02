#!/usr/bin/env bash
set -euo pipefail

# horizontal velcoity file: /home/tao/Prj_GIA_PlateTectonic/Data/plate_motion_presentday/S200_shifted_0 (lon, lat, v_east, v_north)
# note that the first column need shift by 180 degrees to get proper longitude.
# require:
# plot horizontal velocity vectors on a map,
# Using GMT:
# - region extent -130.5/-90.5/20.5/65.5

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/../../.." && pwd)"

vel_file="${repo_dir}/Data/plate_motion_presentday/S200_shifted_0"
out_pdf="${script_dir}/check_plate_motion_pd.pdf"
export GMT_USERDIR="${script_dir}/.gmt"
mag_xyz="${GMT_USERDIR}/plate_motion_magnitude.xyz"
mag_grd="${GMT_USERDIR}/plate_motion_magnitude.nc"
mag_cpt="${GMT_USERDIR}/plate_motion_magnitude.cpt"

region="-130.5/-90.5/20.5/65.5"
projection="Q7i"

# Converts velocity components to GMT vector azimuth/length columns.
# Adjust this scale if the plotted vectors are too long or too short.
vector_scale=0.1
thin_step=2
legend_velocity=4
legend_lon=-128.5
legend_lat=21.5

if [[ ! -f "${vel_file}" ]]; then
    echo "Missing velocity file: ${vel_file}" >&2
    exit 1
fi

mkdir -p "${GMT_USERDIR}"

awk '
    NF >= 4 && $1 !~ /^#/ {
        lon = ($1 + 180.0) % 360.0
        if (lon >= 180.0) {
            lon -= 360.0
        }
        lat = $2
        if (lon >= -130.0 && lon <= -90.0 && lat >= 20.0 && lat <= 65.0) {
            print lon, lat, sqrt($3 * $3 + $4 * $4)
        }
    }
' "${vel_file}" > "${mag_xyz}"

grid_inc="1/1"
mag_range="0/8"

gmt xyz2grd "${mag_xyz}" -G"${mag_grd}" -R"${region}" -I"${grid_inc}"
gmt makecpt -Cwysiwyg -T"${mag_range}" -H > "${mag_cpt}"

gmt begin "${out_pdf%.pdf}" pdf,png
    gmt set MAP_FRAME_TYPE plain FORMAT_GEO_MAP dddF

    gmt coast -R"${region}" -J"${projection}" \
        -Bxa10f5 -Bya10f5 -BWSen+t"Present-day plate motion" \
        -W0.25p,gray40 -N1/0.5p,gray50

    gmt grdimage "${mag_grd}" -R"${region}" -J"${projection}" -C"${mag_cpt}"


    awk -v scale="${vector_scale}" -v thin_step="${thin_step}" '
        NF >= 4 && $1 !~ /^#/ {
            lon = ($1 + 180.0) % 360.0
            if (lon >= 180.0) {
                lon -= 360.0
            }
            lat = $2
            lon_key = sprintf("%.10g", lon)
            lat_key = sprintf("%.10g", lat)
            if (!(lon_key in lon_index)) {
                lon_index[lon_key] = ++n_lon
            }
            if (!(lat_key in lat_index)) {
                lat_index[lat_key] = ++n_lat
            }
            if ((lon_index[lon_key] - 1) % thin_step != 0 || (lat_index[lat_key] - 1) % thin_step != 0) {
                next
            }
            ve = $3
            vn = $4
            speed = sqrt(ve * ve + vn * vn)
            if (speed > 0) {
                azimuth = atan2(ve, vn) * 180.0 / atan2(0, -1)
                if (azimuth < 0) {
                    azimuth += 360.0
                }
                print lon, lat, azimuth, speed * scale
            }
        }
    ' "${vel_file}" | gmt plot -Sv0.18c+e -W0.6p,red -Gred

    printf "%s %s\n" \
        "-129.5 20.8" "-123.0 22.4" \
        | gmt plot -Sr+s -Gwhite -W0.5p,black

    printf "%s %s 0 %s\n" \
        "${legend_lon}" "${legend_lat}" "$(awk -v velocity="${legend_velocity}" -v scale="${vector_scale}" 'BEGIN { print velocity * scale }')" \
        | gmt plot -Sv0.18c+e -W0.6p,red -Gred

    printf "%s %s %s mm/yr\n" \
        "$(awk -v lon="${legend_lon}" 'BEGIN { print lon + 2.5 }')" "${legend_lat}" "${legend_velocity}" \
        | gmt text -F+f10p,Helvetica,black+jML

    gmt colorbar -C"${mag_cpt}" -DJBC+w5i/0.18i+o0/-0.55i+h \
        -Bxaf+l"Velocity magnitude (mm/yr)"
gmt end

echo "Wrote ${out_pdf}"
