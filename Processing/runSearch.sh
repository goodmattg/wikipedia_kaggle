# Compile the script
g++ -Wall DTW_NNSearch.cpp -o DTW_NNSearch

# Run the executable
./DTW_NNSearch \
FrameContainer/DataRows/Data_Rows.txt \
FrameContainer/QueryRows/Query_Rows_0.txt \
550 \
0.01 \
FrameContainer/OutRows/Out_Rows.txt
