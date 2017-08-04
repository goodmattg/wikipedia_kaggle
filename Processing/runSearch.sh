#ioCompile the script
g++ -Wall DTW_NNSearch.cpp -o DTW_NNSearch


END=14
for i in $(seq 1 $END); do

    ./DTW_NNSearch \
    FrameContainer/DataRows/Data_Rows.txt \
    FrameContainer/QueryRows/Query_Rows_${i}.txt \
    550 \
    0.01 \
    FrameContainer/OutRows/Out_Rows_${i}.txt

    echo "Done Processing File ${i}"

done
    # Run the executable
