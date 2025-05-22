# Notes for the C code for the SLE iterations

If previous SVE run is the iteration $k$, now we are preparing surface load file for the next (i.e., $k+1$) SVE run (iteration).



for stage $j$ from 0 ($t_0$) to N ($t_p$):

Now we are doing for time step (stage) $j$.



- topo_presentday: $T_p$
- topo_initial: $T_{j=0}^{k+1}$.
- ice_height_current: $I_j^a$, the original (unmodified) ice height at step j.
- groundIce_mask_current: $G_j$, the ground ice at step j (for last iteration?).
- ocn_current: $O_j^{k+1}$, ocean function to output.
- ocn_initial: $O_0^{k+1}$.
- groundIce_height_current: $I_j^{k+1}$, grounded ice height at step j.
- groundIce_height_initial: $I_0^{k+1}$.
- prevIter_groundIce_height_current: $I_j^{k}$, grounded ice height used in last step.
- prevIter_groundIce_height_initial: $I_0^{k}$.
- prevIter_topo_initial: $T_0^k$, used for $\Delta SL$.
- prevIter_ocn_current: $O_j^k$, ocean function used in the last SVE run.
- prevIter_ocn_initial: $O_{j=0}^k$, ocean function used in the last SVE run.
- prevIter_topo_initial: $T_{j=0}^k$
- $\rho_w$: water density, $\rho_i$: ice density.



In main:

Assign $I_j^{k}$ = $I_j^a (1-O_j^k)$ , also $I_0^k$.

Then in `find_topo_ocean_and_total_load`, for each stage $j$:

1. Subroutine `helper_get_Delta_SL` first gets the $\Delta SL$ at the presetn day, which $=\Delta G_j^k - \Delta R_j^k + c_j^k$, where $c_j^k = -\frac{1}{A_j^k}(\int \frac{\rho_i}{\rho_w} [\Delta {I}_j^k-\frac{\rho_w}{\rho_i}T_0^k(O_j^k-O_0^k)] d\Omega + \int (\Delta G_j^k-\Delta R_j^k)O_j^k d\Omega  )$, where $j=p$. Note 
   1.  here $\Delta I_j^k$ is grounding ice, $=I_j^a (1-O_j^k) - I_0^a (1-O_0^k)$. here also $j=p$.

2. Similar, `helper_get_Delta_SL` get $\Delta SL_j$.
3. Then get updated topography $T_j^{k+1}=T_p^k + \Delta SL^k_p - \Delta SL^k_j$. This means the new initial topography ($T_0^{k+1}$) is also obtained.
4. Then `find_groundice_ocean_mask` using the $T_j^{k+1}$ and $I_j^a$ to check floating ice. So we get new grounding ice $I_j^{k+1}$ and new ocean function $O_j^{k+1}$.
5. Then, `write_ice_and_ocean_func` would write out $I_j^{k+1}-\frac{\rho_w}{\rho_i}T_0^{k+1}(O_j^{k+1}-O_0^{k+1})$ as the surface load (actually relative to the first epoch) for the next SVE iteration, and the $O_j^{k+1}$ for next iteration's ocean functions.

Finally, in `main`, the code would write out the initial topography $T_0^{k+1}$ for next iteration (the $k+2$ -th iteration) to use. 