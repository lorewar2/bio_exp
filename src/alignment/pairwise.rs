
use std::cmp;

pub fn pairwise (seq_x: &Vec<u8>, seq_y: &Vec<u8>, match_score: i32, mismatch_score: i32, gap_open_score: i32, gap_extend_score: i32) -> (Vec<u8>, isize) {
    // variables to save results
    let mut align_vec: Vec<u8> = Vec::new();

    // calculation variables
    // initialize the weight matrix with zeros
    let mut match_matrix: Vec<Vec<isize>> = vec![vec![0; seq_y.len() + 1]; seq_x.len() + 1]; // match or mismatch diagonal edges
    let mut del_matrix: Vec<Vec<isize>> = vec![vec![0; seq_y.len() + 1]; seq_x.len() + 1];  // x deletions right direction edges
    let mut ins_matrix: Vec<Vec<isize>> = vec![vec![0; seq_y.len() + 1]; seq_x.len() + 1]; // x insertion down direction edges
    // initialize the backtrace matrix with ms, ds, and is
    let mut back_matrix: Vec<Vec<char>> = vec![vec!['m'; seq_y.len() + 1]; seq_x.len() + 1];
    back_matrix[0][1] = 'd'; let temp_value = gap_open_score as isize + gap_extend_score as isize; del_matrix[0][1] = temp_value; ins_matrix[0][1] = temp_value; match_matrix[0][1] = temp_value;
    back_matrix[1][0] = 'i'; let temp_value = gap_open_score as isize + gap_extend_score as isize; ins_matrix[1][0] = temp_value; del_matrix[1][0] = temp_value; match_matrix[1][0] = temp_value;
    for i in 2..seq_y.len() + 1 {back_matrix[0][i] = 'd'; let temp_value = del_matrix[0][i - 1] + gap_extend_score as isize; del_matrix[0][i] = temp_value; ins_matrix[0][i] = temp_value; match_matrix[0][i] = temp_value;}
    for i in 2..seq_x.len() + 1 {back_matrix[i][0] = 'i'; let temp_value = ins_matrix[i - 1][0] + gap_extend_score as isize; ins_matrix[i][0] = temp_value; del_matrix[i][0] = temp_value; match_matrix[i][0] = temp_value;}

    // calculations
    // filling out score matrices and back matrix
    for i in 1..seq_x.len() + 1 {
        for j in 1..seq_y.len() + 1 {
            // fill del matrix 
            // get j - 1 score from same matrix with gap extend
            let temp_del_score = del_matrix[i][j - 1] + gap_extend_score as isize;
            // get j - 1 score from match matrix with gap open penalty
            let temp_match_score = match_matrix[i][j - 1] + gap_open_score as isize + gap_extend_score as isize;
            // insert the max
            del_matrix[i][j] = cmp::max(temp_del_score, temp_match_score);

            // fill ins matrix
            // get i - 1 score from the same matrix
            let temp_ins_score = ins_matrix[i - 1][j] + gap_extend_score as isize;
            // get i - 1 score from the match matrix with gap open penalty
            let temp_match_score = match_matrix[i - 1][j] + gap_open_score as isize + gap_extend_score as isize;
            // insert the max
            ins_matrix[i][j] = cmp::max(temp_ins_score, temp_match_score);

            // fill match matrix
            // get the i,j from the insertion matrix
            let temp_ins_score = ins_matrix[i][j];
            // get the i,j from the deletion matrix
            let temp_del_score = del_matrix[i][j];
            // get the match from i-1,j-1 from match matrix with match score or mismatch score
            let temp_match_score;
            if seq_x[i - 1] == seq_y[j - 1] {
                temp_match_score = match_matrix[i - 1][j - 1] + match_score as isize;
                back_matrix[i][j] = 'm';
            }
            else {
                temp_match_score = match_matrix[i - 1][j - 1] + mismatch_score as isize;
                back_matrix[i][j] = 's';
            }
            // insert the max
            match_matrix[i][j] = cmp::max(temp_match_score, cmp::max(temp_ins_score, temp_del_score));
            if (temp_match_score >= temp_ins_score) && (temp_match_score >= temp_del_score) {
                // already allocated
            }
            else if temp_ins_score > temp_del_score {
                back_matrix[i][j] = 'i';
            }
            else {
                back_matrix[i][j] = 'd';
            }
        }
    }
    // print stuff
    println!("simple ins matrix");
    for i in 0..seq_x.len() + 1 {
        for j in 0..seq_y.len() + 1 {
            print!("{:>3} ", ins_matrix[i][j]);
        }
        println!("");
    }
    println!("simple match matrix");
    for i in 0..seq_x.len() + 1 {
        for j in 0..seq_y.len() + 1 {
            print!("{:>3} ", match_matrix[i][j]);
        }
        println!("");
    }
    println!("simple del matrix");
    for i in 0..seq_x.len() + 1 {
        for j in 0..seq_y.len() + 1 {
            print!("{:>3} ", del_matrix[i][j]);
        }
        println!("");
    }
    println!("simple back matrix");
    for i in 0..seq_x.len() + 1 {
        for j in 0..seq_y.len() + 1 {
            print!("{:>3} ", back_matrix[i][j]);
        }
        println!("");
    }
    // back tracing using back matrix and filling out align_vec
    let mut i = seq_x.len();
    let mut j = seq_y.len();
    let score = match_matrix[i][j];
    let mut break_on_next = false;
    loop {
        match back_matrix[i][j] {
            'i' => {
                i = i - 1;
                align_vec.push('i' as u8);
            },
            'm' => {
                i = i - 1;
                j = j - 1;
                align_vec.push('m' as u8);
            },
            's' => {
                i = i - 1;
                j = j - 1;
                align_vec.push('s' as u8);
            }
            'd' => {
                j = j - 1;
                align_vec.push('d' as u8);
            },
            _ => (),
        }
        if break_on_next {
            break;
        }
        if i == 0 && j == 0 {
            break_on_next = true;
        }
    }
    (align_vec.into_iter().rev().collect(), score)
}