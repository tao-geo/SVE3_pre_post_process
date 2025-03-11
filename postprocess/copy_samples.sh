# copy some samples from combined_files_case_test to combined_files_case_test_samples

step_list=(16  768  856  976)  # first output, 26 ka, 15 ka, 0 ka

for step in ${step_list[@]}; do
    cp combined_files_case_test/*.${step}.*     combined_files_case_test_samples/
    cp combined_files_case_test/*.${step}       combined_files_case_test_samples/
done

# rm -r ./combined_files_case_test