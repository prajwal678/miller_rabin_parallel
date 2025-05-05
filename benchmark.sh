#!/bin/bash
BIT_SIZES=(512 1024 2048 4096)
ITERATIONS=50
TEST_REPETITIONS=3

if ! command -v bc &> /dev/null; then
    echo "error: bc is not there."
    exit 1
fi
if ! command -v xxd &> /dev/null; then
    echo "error: xxd is not there."
    exit 1
fi

echo "Compiling code right now
make clean && make
if [ $? -ne 0 ]; then
    echo "compilation failed, exiting."
    exit 1
fi

SEQ_EXE="bin/miller_rabin_seq"
PAR_EXE="bin/miller_rabin_par"

if [ ! -f "$SEQ_EXE" ]; then
    echo "error: sequential executable not found at $SEQ_EXE"
    exit 1
fi

mkdir -p results

generate_random_number() { # ze for the testing
    local bits=$1
    local bytes=$(( (bits + 7) / 8 ))
    
    # byte to hex
    local random_hex=$(dd if=/dev/urandom bs=1 count=$bytes 2>/dev/null | xxd -p -c $bytes)
    
    # check if MSB is set for correct bit length
    local first_byte_hex="0x${random_hex:0:2}"
    local first_byte_dec=$((first_byte_hex))
    local msb_in_first_byte_pos=$(( (bits - 1) % 8 ))
    local msb_mask=$(( 1 << msb_in_first_byte_pos ))
    local first_byte_dec_msb=$(( first_byte_dec | msb_mask ))
    local first_byte_hex_msb=$(printf "%02x" $first_byte_dec_msb)
    random_hex="${first_byte_hex_msb}${random_hex:2}"

    # check if LSB is set (make it odd)
    local last_byte_hex="0x${random_hex: -2}"
    local last_byte_dec=$((last_byte_hex))
    local last_byte_dec_lsb=$(( last_byte_dec | 1 ))
    local last_byte_hex_lsb=$(printf "%02x" $last_byte_dec_lsb)
    if [ ${#random_hex} -gt 2 ]; then
        random_hex="${random_hex:0:-2}${last_byte_hex_lsb}"
    else
        random_hex="$last_byte_hex_lsb"
    fi

    # convert final hex to decimal using bc, removing any line breaks ( yes these shiz are LLM-ed)
    local random_dec=$(echo "ibase=16; ${random_hex^^}" | bc | tr -d '\\ 
')

    if [ -z "$random_dec" ]; then
        echo "error: bc failed to convert hex to dec" >&2
        echo "Hex value was: $random_hex" >&2
        return 1 
    fi

    echo "$random_dec"
}

run_test() {
    local bits=$1
    local iterations=$2
    local repetitions=$3
    
    echo "tests for $bits-bit numbers with $iterations iterations:->"
    
    local test_number=$(generate_random_number $bits)
    if [ $? -ne 0 ] || [ -z "$test_number" ]; then
        echo "error: generating random number for $bits bits. skipping test."
        printf "| %-8s | %-10s | %-17s | %-15s | %-7s |\n" "$bits" "$iterations" "Error" "Error" "N/A"
 >> results/summary.md
        return
    fi
    echo "  testing number (first 50 digits): ${test_number:0:50}..."
    
    local seq_total=0
    local par_total=0
    local speedup_total=0
    local valid_repetitions=0
    local valid_speedup_calcs=0
    local par_available=0
    if [ -f "$PAR_EXE" ]; then
        par_available=1
    fi
    
    for ((i=1; i<=repetitions; i++)); do
        echo "    Repetition $i of $repetitions"
        local seq_time_ok=0
        local par_time_ok=0
        local current_seq_time=""
        local current_par_time=""
        
        echo "      running sequential test..."
        SEQ_OUTPUT=$($SEQ_EXE "$test_number" $iterations)
        current_seq_time=$(echo "$SEQ_OUTPUT" | grep -o 'Total execution time: [0-9.]* ms' | awk '{print $4}')
        # Validate time format (basic check for digits and dot), having to do this cuz GOOGLE COLAB HAS CLOCK SKEW FOR SM RSN AUGHHH
        if [[ "$current_seq_time" =~ ^[0-9.]+$ ]]; then
            echo "      Sequential time: $current_seq_time ms"
            seq_total=$(echo "$seq_total + $current_seq_time" | bc)
            seq_time_ok=1
        else
            echo "      Error extracting valid time from sequential output."
            echo "      Output: $SEQ_OUTPUT"
            current_seq_time="Error" # Mark as error
        fi
        
        if [ $par_available -eq 1 ]; then
            echo "      running parallel test..."
            PAR_OUTPUT=$($PAR_EXE "$test_number" $iterations)
            current_par_time=$(echo "$PAR_OUTPUT" | grep -o 'Total execution time: [0-9.]* ms' | awk '{print $4}')
            # Validate time format
            if [[ "$current_par_time" =~ ^[0-9.]+$ ]]; then
                echo "      Parallel time: $current_par_time ms"
                par_total=$(echo "$par_total + $current_par_time" | bc)
                par_time_ok=1
                
                if [ $seq_time_ok -eq 1 ] && [ $(echo "$current_seq_time > 0" | bc) -eq 1 ] && [ $(echo "$current_par_time > 0" | bc) -eq 1 ]; then
                    # higher scale for division first, then format
                    current_speedup=$(echo "scale=4; $current_seq_time / $current_par_time" | bc)
                    current_speedup_fmt=$(printf "%.2f" $current_speedup)
                    echo "      Speedup: $current_speedup_fmt x"
                    speedup_total=$(echo "$speedup_total + $current_speedup" | bc) # Use unformatted for sum
                    valid_speedup_calcs=$((valid_speedup_calcs + 1))
                elif [ $seq_time_ok -eq 1 ]; then
                    echo "      Speedup: N/A (Seq or Par time was zero or invalid)" 
                fi
            else
                echo "      error extracting valid time from parallel output."
                echo "      Output: $PAR_OUTPUT"
                current_par_time="Error" # Mark as error
            fi
        else
            echo "      parallel executable not found. skipping parallel test."
            par_time_ok=1 # Consider it "ok" for averaging if not available
        fi
        
        # count repetition only if both seq and par (if avail) extracted valid times
        if [ $seq_time_ok -eq 1 ] && ([ $par_available -eq 0 ] || [ $par_time_ok -eq 1 ]); then
             valid_repetitions=$((valid_repetitions + 1))
        else
            echo "    Skipping repetition $i from average due to time extraction error."
        fi
    done
    
    echo "Debug: Valid repetitions for $bits bits: $valid_repetitions"
    if [ $valid_repetitions -gt 0 ]; then
        seq_avg=$(echo "scale=2; $seq_total / $valid_repetitions" | bc)
        [[ "$seq_avg" == .* ]] && seq_avg="0$seq_avg"

        if [ $par_available -eq 1 ]; then
            par_avg=$(echo "scale=2; $par_total / $valid_repetitions" | bc)
            [[ "$par_avg" == .* ]] && par_avg="0$par_avg"
            if [ $valid_speedup_calcs -gt 0 ]; then
                 speedup_avg_raw=$(echo "scale=4; $speedup_total / $valid_speedup_calcs" | bc)
                 speedup_avg=$(printf "%.2f" $speedup_avg_raw)
            else 
                 speedup_avg="N/A"
            fi
        else
            par_avg="N/A"
            speedup_avg="N/A"
        fi
        
        echo "  Average results for $bits-bit numbers ($valid_repetitions valid repetitions):"
        echo "    Sequential: $seq_avg ms"
        echo "    Parallel: $par_avg ms"
        echo "    Speedup: $speedup_avg x"
        
        row_output=$(printf "| %-8s | %-10s | %-17s | %-15s | %-7s |" "$bits" "$iterations" "$seq_avg" "$par_avg" "$speedup_avg")
        echo "Debug: Appending row (using echo): $row_output"
        echo "$row_output" >> results/summary.md
        if [ $? -ne 0 ]; then
            echo "Debug: error appending to results/summary.md using echo"
        fi
    else
        echo "  No valid repetitions completed for $bits-bit numbers. Cannot calculate averages."
        error_row_output=$(printf "| %-8s | %-10s | %-17s | %-15s | %-7s |" "$bits" "$iterations" "Error" "Error" "N/A")
        echo "Debug: Appending error row (using echo): $error_row_output"
        echo "$error_row_output" >> results/summary.md
        if [ $? -ne 0 ]; then
            echo "Debug: Error appending error row to results/summary.md using echo"
        fi
    fi
}

echo "# Miller-Rabin Primality Test Performance Comparison" > results/summary.md
echo "" >> results/summary.md
echo "Configuration:" >> results/summary.md
echo "- Miller-Rabin iterations: $ITERATIONS" >> results/summary.md
echo "- Test repetitions per bit size: $TEST_REPETITIONS" >> results/summary.md
echo "- Random number generated per bit size using /dev/urandom" >> results/summary.md
echo "" >> results/summary.md
printf "| %-8s | %-10s | %-17s | %-15s | %-7s |\n" "Bit Size" "Iterations" "Sequential (ms)" "Parallel (ms)" "Speedup"

printf "|%s|%s|%s|%s|%s|\n" "----------" "------------" "-----------------" "---------------" "---------|" >> results/summary.md

for bits in "${BIT_SIZES[@]}"; do
    run_test $bits $ITERATIONS $TEST_REPETITIONS
done

echo ""
echo "Benchmark completed. Results saved to results/summary.md"
echo ""
echo "Summary:"
cat results/summary.md