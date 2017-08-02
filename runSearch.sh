# Compile the script
g++ -Wall UCR-DTW/SourceCode/UCR_DTW_MOD.cpp  -o UCR-DTW/Executable/UCR_DTW_MOD

# Run the executable
./UCR-DTW/Executable/UCR_DTW_MOD Processing/Data_Rows.txt Processing/Data_Rows.txt 550 0.01 Processing/Out_Rows.txt
