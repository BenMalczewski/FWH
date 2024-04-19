# FWH
This github contains all of the required files to install a compatable version of Julia along with the required packages to run the FWH code. This code is confirmed to function correctly on native Linux workstations and Ubuntu WSL. Using this code on an HPC cluster has yet to be successful, due to incompatabilities with MPI. If installing locally, it is recommended to use intel's oneapi.
https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html?cid=sem&source=sa360&campid=2024_ao_cbu_us_gmo_gmocrbu_awa_text-link_brand_exact_cd_ai-oneapi_3500186468_google_b2b_is_non-pbm_intel&ad_group=AI_Brand-oneAPI_Oneapi_Exact&intel_term=intel+oneapi&sa360id=43700079829610580&gad_source=1&gclid=CjwKCAjwrIixBhBbEiwACEqDJdEmcmQEqb8d-3HOrLrSfb-dDMlED_eHXNWGendQlx-kDB7M_Km0yBoChBEQAvD_BwE&gclsrc=aw.ds#gs.876q61

Download the installJulia.zip folder. Extract the contents into your $HOME directory. 
```bash
cd installJulia
chmod +x *
./installJulia
```

