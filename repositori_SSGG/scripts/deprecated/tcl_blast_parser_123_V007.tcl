#!/usr/bin/tclsh

##############################################################
#                                                            #
#               Tcl/Tk NCBI BLAST PARSER                     #
#              COPYRIGHT, Alexander Kozik                    #
#                                                            #
##############################################################
#                                                            #
# FROM BLAST OUTPUT PARSER EXTRACTS:                         #
#                                                            #
#    1. Query ID                                             #
#    2. First best hit ID                                    #
#    3. Normalized Expectation: (-log(exp_value))            #
#       if "exp" less than 1e-100 value 1e-100 assigned      #
#    4. Identity(%)                                          #
#    5. Hit number                                           #
#           Data 1 - 5 are written into "best_hit" file      #
#                                                            #
#    6. Total number of hits                                 #
#       Data 6 are written into separate "blast_stat" file   #
#                                                            #
#    7. Total list of all hits (primary and alternative      #
#       alignments).                                         #
#       Data 7 are written into separate "all_hits" file     #
#                                                            #
#    8. "Tho best hits" are written into "two_hits" file     #
#       as well as info about their expectation difference   #
#                                                            #
#    9. List of all "genes" in BLAST output is written       #
#       to "ID_list" file (we need it to create              #
#       true "Matrix" file for GenomePixelizer or            #
#       PhyloGrapher.                                        #
#                                                            #
#   10. Matrix file (matrix). Where third column is identity #
#       value (in the range from 0 to 1) between pairs       # 
#       of genes (nodes), fourth column is expectation value #
#       and fifth column is overlap size.                    #
#                                                            #
#   11. Graph groups (in the directory "subgroups"):         #
#       Sets of genes (nodes) connected to each other        #
#       based on identity cutoff value using DFS             #
#       (Depth First Search) algorithm by analysis of        #
#       the adjacency list file ("adj_list" extension)       #
#                                                            #
#   12. Group Degree Info file ("group_degree_info"          #
#       extension) with info about every gene (node)         #
#       belonging to distinct group, adjacency number,       #
#       group complexity and group number (ID)               #
#       based on graph analysis.                             #
#                                                            #
##############################################################
#                                                            #
#  PROGRAM USAGE:                                            #
#                                                            #
# tcl_blast_parser_123.tcl [plus six arguments]              #
#                                                            #
#  ARGUMENTS:                                                #
#                                                            #
# 1. input (BLAST results)                                   #
# 2. output (output file name)                               #
# 3. expectation cutoff value                                #
#    (hits only with this "exp" value or better will be used #
#     to generate matrix file and in full graph analysis)    #
# 4. identity cutoff value                                   #
#    (hits only with this identity or better will be used    #
#     to generate matrix file and in full graph analysis)    #
# 5. overlap cutoff value                                    #
#    (hits only with this overlap or better will be used     #
#     to generate matrix file and in full graph analysis)    #
# 6. MATRIX/GRAPH option to choose                           #
#    where to stop BLAST file processing:                    #
#       A. on matrix file                                    #
#              or on                                         #
#       B. full graph analysis                               #
#    (graph analysis may work really long time with          #
#    large set of genes (more than 10,000)                   #
#                                                            #
#  FOR EXAMPLE:                                              #
#                                                            ########
#                                                                   #
# tcl_blast_parser_123.tcl my_blast.out my_results 20 40 100 GRAPH  #
#                                                                   #
# In this case program will read and analyse "my_blast.out"  ########
# file and generate eight output files plus one additional   #
# directory "subgroups" with list of groups of genes based   #
# on clustering analysis.                                    ########
#                                                                   #
# tcl_blast_parser_123.tcl my_blast.out my_results 20 40 100 MATRIX #
#                                                                   #
# In this case program will read and analyse "my_blast.out"  ########
# file and generate six output files skipping detailed       #
# graph analysis.                                            #
#                                                            #
##############################################################


    #### FUNCTIONS ####

    proc Blast_Parsing {argv} {

	set f_in [open [lindex $argv 0] "r"]
	set f_out1 [open [lindex $argv 1].best_hit "w"]
	set f_out2 [open [lindex $argv 1].blast_stat "w"]
	set f_out3 [open [lindex $argv 1].all_hits "w"]
	set f_out4 [open [lindex $argv 1].two_hits "w"]
	set f_out5 [open [lindex $argv 1].id_list "w"]

	### "q" IS A COUNTER ###
	### SET "q" COUNTER TO 1 ###
	set q  1
	### BLANK LINE COUNTER ###
	set blank_line_counter 0
	### TRICK WITH SUBJECT END ###
	set query_end_find "NOT_YET"
	set subject_end_find "NOT_YET"
	### THIRD ARGUMENT TO DISPLAY ALT ALIGNMENT ###
	set display_alt_align "TRUE"
	### CUTOFF VALUES ###
	# set exp_cut 20 ; # RECOMMENDED DEFAULT
	# set idn_cut 40 ; # RECOMMENDED DEFAULT
	# set ovr_cut 100 ; # RECOMMENDED DEFAULT
	set exp_cut [lindex $argv 2]
	set idn_cut [lindex $argv 3]
	set ovr_cut [lindex $argv 4]
	set graph_analysis [lindex $argv 5]
	puts "$exp_cut Expect Cutoff"
	puts "$idn_cut Identity Cutoff"
	puts "$ovr_cut Overlap Cutoff"
	puts "Stop analysis at: $graph_analysis"
	puts " Starting BLAST Parsing"
	puts "======================="
	after 1000

	while {[gets $f_in t] >= 0} {

	    if {$t == ""} {
		incr blank_line_counter
		# puts "$blank_line_counter\BLANK_LINE"
	    }

	    set query_match [string range $t 0 6]
	    if {$query_match == "Query= "} {
		### UNSET "query_id" FOR DEBUGGING PURPOSE ###
		set query_id "WHATEVER"
		### UNSET "subject_found" FOR DEBUGGING PURPOSE ###
		set subject_found "NO"
		### UNSET "e" (EXPECT) FOR DEBUGGING PURPOSE ###
		set e "WHATEVER"
		### UNSET QUERY AND SUBJECT LENGTH ###
		set query_start_find "NOT_YET"
		set subject_start_find "NOT_YET"
		set query_end_find "NOT_YET"
		set subject_end_find "NOT_YET"
		### COUNT NUMBER OF HITS ###
		### SET "hit_counter" TO "0" ###
		set hit_counter 0
		set two_hits_found "FALSE"
		regsub {^Query= } $t "" qr
		regsub {\(} $qr "" qr
		regsub {\)} $qr "" qr
		regsub {\,} $qr "" qr
		set query_id $qr
		regsub {^gi\|} $query_id "" query_id
		regsub {\|.*} $query_id "" query_id
		regsub { .*} $query_id "" query_id
		puts "$query_id\t$q"
		puts -nonewline $f_out2 "$q\t$query_id\t"
		if {$q == 1} {
		    ### LIST OF GENES IN BLAST OUTPUT INITIALIZATION ###
		    set gene_list [lappend gene_list $query_id]
		} else {
		    set already_done [lsearch -exact $gene_list $query_id]
		    if {$already_done < 0} {
			set gene_list [lappend gene_list $query_id]
		    }
		}
		incr q
	    }

	    set subject_match [string range $t 0 0]
	    if {$subject_match == ">"} {
		### UNSET ALT ALIGNMENT ###
		set alt_align_count 0
		regsub {^>} $t "" subj
		regsub {^gi\|} $subj "" subj
		regsub {\|.*} $subj "" subj
		regsub { .*} $subj "" subj
		set already_done [lsearch -exact $gene_list $subj]
		if {$already_done < 0} {
		    set gene_list [lappend gene_list $subj]
		}
		incr hit_counter
		set subject_found "YES"
	    }

	    set score_match [string range $t 0 8]
	    if {$score_match == " Score = "} {
		### UNSET QUERY AND SUBJECT LENGTH ###
		set query_start_find "NOT_YET"
		set subject_start_find "NOT_YET"
		set query_end_find "NOT_YET"
		set subject_end_find "NOT_YET"
		### ALT ALIGN COUNTER ###
		incr alt_align_count
		if {$alt_align_count > 1} {
		    puts "$query_id\t$subj\talt_alignment"
		}
		### TBLASTX CASE ###
		regsub {.*Expect\([0-9]\) = } $t "" t
		regsub {.*Expect\([0-9][0-9]\) = } $t "" t
		### REGULAR CASE ###
		regsub {.*Expect = } $t "" t
		if {[string range $t 0 0] == "e"} {
		    regsub {e} $t "1e" t
		}
		if {$subject_found == "YES"} {
		    set e $t
		    if {$e == 0.0} {
			set e 1e-254
		    }
		    set e [expr -(log10($e))]
		    if {$e == "inf"} {
			set e 100
		    }
		    if {$e > 100} {
			set e 100
		    }
		    set e [expr ceil($e)]
		    set e [expr int($e)]
		    if {$hit_counter == 1 && $alt_align_count == 1} {
			set first_hit_id $subj
			set first_hit_exp $e
		    }
		    if {$hit_counter == 2 && $alt_align_count == 1} {
			set second_hit_id $subj
			set second_hit_exp $e
			set two_hits_found "TRUE"
		    }
		}
	    }

	    set identity_match [string range $t 0 13]
	    if {$identity_match == " Identities = "} {
		regsub {\, Gaps =.*} $t "" idn
		regsub {\, Positives =.*} $idn "" idn
		regsub {.*\(} $idn "" idn
		regsub {\%\)} $idn "" idn
		regsub { Identities = } $t "" frc
		regsub { .*} $frc "" frc
		regsub {\/} $frc "\t" frc

		if {$hit_counter == 1 && $alt_align_count == 1} {
		    set first_hit_idn $idn
		    set first_hit_frc [lindex [split $frc] 1]
		}
		if {$hit_counter == 2 && $alt_align_count == 1} {
		    set second_hit_idn $idn
		    set second_hit_frc [lindex [split $frc] 1]
		}

		if {$subject_found == "YES"} {
		    # set subject_found "NO"
		    if {$alt_align_count == 1} {
			puts -nonewline $f_out3 "$query_id\t$subj\t$e\t$idn\t$frc\t$hit_counter\t$alt_align_count\tPRM\t"
		    }
		    if {$alt_align_count > 1 && $display_alt_align == "TRUE"} {
			puts -nonewline $f_out3 "$query_id\t$subj\t$e\t$idn\t$frc\t$hit_counter\t$alt_align_count\tALT\t"
		    }
		    if {$hit_counter == 1 && $alt_align_count == 1} {
			puts $f_out1 "$query_id\t$subj\t$e\t$idn\t$frc\t$hit_counter"
		    }
		}
	    }

	    set frame_match [string range $t 0 8]
	    if {$frame_match == " Frame = "} {
		regsub { Frame = } $t "" fr
		regsub -all { } $fr "" fr
		puts -nonewline $f_out3 "$fr\t"
	    }

	    set strand_match [string range $t 0 9]
	    if {$strand_match == " Strand = "} {
		regsub { Strand = } $t "" str
		regsub -all {Plus} $str "+" str
		regsub -all {Minus} $str "-" str
		regsub -all { } $str "" str
		puts -nonewline $f_out3 "$str\t"
	    }

	    set query_start [string range $t 0 5]
	    if {$query_start == "Query:" && $query_start_find == "NOT_YET"} {
		regsub {^Query\: } $t "" qs
		regsub { .*} $qs "" qs
		# puts $qs
		set query_start_find "DONE"
	    }

	    set subject_start [string range $t 0 5]
	    if {$subject_start == "Sbjct:" && $subject_start_find == "NOT_YET"} {
		regsub {^Sbjct\: } $t "" ss
		regsub { .*} $ss "" ss
		# puts $ss
		set subject_start_find "DONE"
	    }

	    set query_end [string range $t 0 5]
	    if {$query_end == "Query:"} {
		regsub {.* } $t "" qe
		# puts $qe
		set query_end_find "DONE"
		set blank_line_counter 0
	    }

	    set subject_end [string range $t 0 5]
	    if {$subject_end == "Sbjct:"} {
		regsub {.* } $t "" se
		# puts $se
		set subject_end_find "DONE"
	    }

	    if {$blank_line_counter == 2 && $query_end_find == "DONE" && $subject_end_find == "DONE"} {
		# puts "$qe\t$se\tthis is the END!\tblank_line_counter"
		puts $f_out3 "$qs\t$qe\t$ss\t$se"
		set query_end_find "WRITTEN"
		set subject_end_find "WRITTEN"
	    }

	    ### STAND ALONE BLAST ###
	    set no_hits_match [string range $t 0 26]
	    if {$no_hits_match == " ***** No hits found ******"} {
		# set t "No hits found"
		puts $f_out3 "$query_id\tno_hits_found\t0\t0"
		puts $f_out1 "$query_id\tno_hits_found\t0\t0"
	    }

	    ### NCBI WEB BLAST ###
	    set no_hits_match [string range $t 0 30]
	    if {$no_hits_match == "No significant similarity found"} {
		# set t "No hits found"
		puts $f_out3 "$query_id\tno_hits_found\t0\t0"
		puts $f_out1 "$query_id\tno_hits_found\t0\t0"
	    }

	    set end_of_blast [string range $t 0 32]
	    if {$end_of_blast  == "  Number of sequences in database" && $hit_counter != 1} {
		puts $f_out2 "$hit_counter"
	    }
	    if {$end_of_blast  == "  Number of sequences in database" && $hit_counter == 1 && $e >= 20} {
		puts $f_out2 "$hit_counter\tSTRONG_SINGLE_HIT"
	    }
	    if {$end_of_blast  == "  Number of sequences in database" && $hit_counter == 1 && $e < 20} {
		puts $f_out2 "$hit_counter\tWEAK_SINGLE_HIT"
	    }
	    if {$end_of_blast  == "  Number of sequences in database" && $two_hits_found == "TRUE"} {
		set exp_diff [expr $first_hit_exp - $second_hit_exp]
		set idn_diff [expr ($first_hit_idn*1.00) / ($second_hit_idn*1.00)]
		set frc_diff [expr $first_hit_frc - $second_hit_frc]
		set diff_quality "BAD_DIFF"

		##########################################################################################
		###                     PARAMETERS FOR "TWO_HITS" FILE                                 ###
		###          MODIFY CUTOFF VALUES FOR YOUR OWN PARTICULAR PROJECT NEEDS                ###
		###                                                                                    ###
		### DEFAULT SETTINGS:                                                                  ###
		###                                                                                    ###
		### "exp_diff_cutoff" - SECOND HIT HAS EXPECTATION DIFF 1e-20 OR BETTER                ###
		### "idn_diff_cutoff" - SECOND HIT HAS IDENTITY TWO TIMES LESS THAN FIRST HIT          ###
		### "frc_diff_cutoff" - SECOND HIT HAS OVERLAP (ALIGNMENT) LENGTH LESS THAN FIRST HIT  ###
		###                                                                                    ###
		### IN THIS CASE SECOND HIT QUALIFIED AS A HIT WITH "GOOD_DIFF" IN "TWO_HITS FILE      ###
		##########################################################################################

		set exp_diff_cutoff 20 ; ## <- MODIFY EXPECTATION VALUE
		set idn_diff_cutoff 2  ; ## <- MODIFY IDENTITY VALUE
		set frc_diff_cutoff 0  ; ## <- MODIFY OVERLAP (ALIGNMENT) LENGTH VALUE

		#### END OF MODIFICATIONS ####

		if {$exp_diff >= $exp_diff_cutoff && $idn_diff >= $idn_diff_cutoff && $frc_diff >= $frc_diff_cutoff} {
		    set diff_quality "GOOD_DIFF"
		}
		puts $f_out4 "$query_id\t$first_hit_id\t$first_hit_exp\t$first_hit_idn\t$first_hit_frc"
		puts $f_out4 "$query_id\t$second_hit_id\t$second_hit_exp\t$second_hit_idn\t$second_hit_frc\t$exp_diff\t$idn_diff\t$frc_diff\t$diff_quality"
	    }

	}

	# puts $gene_list
	set gene_list [lsort $gene_list]
	# puts $gene_list

	set k 0

	puts ""
	foreach gene $gene_list {
	    puts $f_out5 $gene
	    set k_mod [expr fmod($k,10)]
	    if {$k_mod == 0} {
		puts -nonewline "."
	    }
	    incr k
	}
	puts ""

	close $f_in
	close $f_out1
	close $f_out2
	close $f_out3
	close $f_out4
	close $f_out5

	Redundant_Matrix $argv $exp_cut $idn_cut $ovr_cut $graph_analysis

    }

    proc Redundant_Matrix {argv exp_cut idn_cut ovr_cut graph_analysis} {

	global matrix_array

	set f_in2 [open [lindex $argv 1].all_hits "r"]
	set f_out6 [open [lindex $argv 1].matrix "w"]

	# puts [lindex $argv 1].all_hits
	puts "Starting Matrix Extraction"
	puts "=========================="
	after 1000
	puts ""
	puts "Reading \"All Hits\" file"

	set l 0

	while {[gets $f_in2 current_line] >= 0}  {

	    set current_line [split $current_line]

	    set id_A [lindex $current_line 0]
	    set id_B [lindex $current_line 1]
	    set expt [lindex $current_line 2]
	    set idnt [lindex $current_line 3]
	    set ovrl [lindex $current_line 5]
	    set algn [lindex $current_line 8]

	    # puts "$id_A\t$id_B\t$expt\t$idnt\t$ovrl\t$algn"

	    if {$id_A != $id_B && $algn == "PRM" && $expt >= $exp_cut && $idnt >= $idn_cut && $ovrl >= $ovr_cut} {

		set matrix_array($id_A,$id_B) "$expt $idnt $ovrl"
		set matrix_list [lappend matrix_list $matrix_array($id_A,$id_B)]
		set pairs_list [lappend pairs_list "$id_A $id_B"]
		set l_mod [expr fmod($l,100)]
		if {$l_mod == 0} {
		    puts -nonewline "."
		}
		incr l
	    }
	}

	puts ""
	puts "Processing Matrix"

	foreach pair $pairs_list {
	    set pair [split $pair]
	    set id_A [lindex $pair 0]
	    set id_B [lindex $pair 1]
	    set query1 [info exists matrix_array($id_A,$id_B)]
	    # puts $query1
	    set query2 [info exists matrix_array($id_B,$id_A)]
	    # puts $query2
	    if {$query1 == 1} {
		set matrix_values1 [split $matrix_array($id_A,$id_B)]
		set exp_value1 [lindex $matrix_values1 0]
		set idn_value1 [lindex $matrix_values1 1]
		set ovr_value1 [lindex $matrix_values1 2]
	    }
	    if {$query2 == 1} {
		set matrix_values2 [split $matrix_array($id_B,$id_A)]
		set exp_value2 [lindex $matrix_values2 0]
		set idn_value2 [lindex $matrix_values2 1]
		set ovr_value2 [lindex $matrix_values2 2]
	    }
	    # puts "$id_A\t$id_B\t$exp_value1\t$idn_value1\t$ovr_value1"
	    # puts "$id_B\t$id_A\t$exp_value2\t$idn_value2\t$ovr_value2"
	    if {$query1 == 1 && $query2 == 0} {
		set idn_value1 [expr $idn_value1/100.00]
		puts $f_out6 "$id_A\t$id_B\t$idn_value1\t$exp_value1\t$ovr_value1"
		unset matrix_array($id_A,$id_B)
	    }
	    if {$query2 == 1 && $query1 == 0} {
		set idn_value2 [expr $idn_value2/100.00]
		puts $f_out6 "$id_A\t$id_B\t$idn_value2\t$exp_value2\t$ovr_value2\tERROR"
		puts "Something is wrong..."
		unset matrix_array($id_B,$id_A)
	    }
	    if {$query1 == 1 && $query2 == 1} {
		if {$idn_value1 >= $idn_value2} {
		    set idn_value1 [expr $idn_value1/100.00]
		    set idn_value2 [expr $idn_value2/100.00]
		    puts $f_out6 "$id_A\t$id_B\t$idn_value1\t$exp_value1\t$ovr_value1"
		}
		if {$idn_value1 < $idn_value2} {
		    set idn_value1 [expr $idn_value1/100.00]
		    set idn_value2 [expr $idn_value2/100.00]
		    puts $f_out6 "$id_A\t$id_B\t$idn_value2\t$exp_value2\t$ovr_value2"
		}
		unset matrix_array($id_A,$id_B)
		unset matrix_array($id_B,$id_A)
	    }
	    set l_mod [expr fmod($l,100)]
	    if {$l_mod == 0} {
		puts -nonewline "."
	    }
	    incr l
	}

	puts ""
	puts "=========================="
	puts "| Matrix Extraction Done |"
	puts "=========================="
	close $f_in2
	close $f_out6

	if {$graph_analysis == "GRAPH"} {
	    Graph_Extraction $idn_cut $argv
	} else {
	    exit
	}
    }

    proc Graph_Extraction {idn_cut argv} {

	puts ""
	puts "Group Extraction Begin"
	puts ""

	after 1000

	global matrix_array

	set f_in3 [open [lindex $argv 1].id_list "r"]
	set f_in4 [open [lindex $argv 1].matrix "r"]

	### CREATE LIST WITH NODE IDs ###
	set k 1
	while {[gets $f_in3 current_line] >= 0}  {
	    set id_list [lappend id_list $current_line]
	    puts $k\t$current_line
	    incr k
	}

	### CREATE MATRIX ARRAY ###
	puts "Reading Matrix File"
	set l 1
	while {[gets $f_in4 current_line] >= 0}  {
	    set current_line [split $current_line]
	    set id_A [lindex $current_line 0]
	    set id_B [lindex $current_line 1]
	    set idnt [lindex $current_line 2]
	    set matrix_array($id_A,$id_B) $idnt
	    set pair_list [lappend pair_list "$id_A $id_B"]
	    set l_mod [expr fmod($l,100)]
	    if {$l_mod == 0} {
		puts -nonewline "."
	    }
	    incr l
	}

	puts ""

	Adj_List_Extraction $idn_cut $argv $id_list $pair_list

	puts ""
	puts "End of Round 1 (Ajacency List)"
	close $f_in3
	close $f_in4
	puts "Starting DFS Procedure"
	after 1000
	DFS_Algorithm $argv $id_list
    }

    proc Adj_List_Extraction {idn_cut argv id_list pair_list} {

	global matrix_array
	global adj_array

	set f_out7 [open [lindex $argv 1].adj_list "w"]

	set idn_cut [expr $idn_cut/100]

	foreach id $id_list {
	    puts $id
	    if {[info exists adj_list] == 1} {
		unset adj_list
	    }
	    set adj_list [lappend adj_list $id]
	    foreach pair $pair_list {
		set pair [split $pair]
		set a [lindex $pair 0]
		set b [lindex $pair 1]
		set value $matrix_array($a,$b)
		# puts $a\t$b\t$value
		if {$id == $a && $value >= $idn_cut} {
		    set adj_list [lappend adj_list $b]
		}
		if {$id == $b && $value >= $idn_cut} {
		    set adj_list [lappend adj_list $a]
		}
		if {$id != $a && $id != $b} {
		    continue
		}
	    }
	    set adj_array($id) $adj_list
	    puts $f_out7 $adj_list
	}
	close $f_out7
    }

    proc DFS_Algorithm {argv id_list} {

	global adj_array
	global used_nodes_list
	global group_members_counter
	global round_2_counter
	global retarded_counter
	global matrix_array
	global sub_group_list
	global group_array
	global group_stat_array
	global list_of_group_array

	set round_2_counter 1

	### CREATE DIRECTORY WITH SUBGROUPS ###
	file delete -force "subgroups"
	file mkdir "subgroups"

	set m 0
	foreach id $id_list {
	    if {$m == 0} {
		set used_nodes_list [lappend used_nodes_list ""]
	    }
	    set already_done [lsearch -exact $used_nodes_list $id]
	    if {$already_done < 0} {
		set sub_group_list ""
		set used_nodes_list [lappend used_nodes_list "$id"]
		incr m
		set sub_group "subgroup_$m"
		set sub_matrx "submatrix_$m"
		set sub_group_open [open "subgroups\/$sub_group" "w"]
		set sub_matrx_open [open "subgroups\/$sub_matrx" "w"]
		set group_members_counter 1
		set group_array($id) "$id\t$m\t*****"
		set list_of_group_array [lappend list_of_group_array $group_array($id)]
		puts $sub_group_open $id
		set sub_group_list [lappend sub_group_list $id]
		set retarded_counter 1
		DFS_Procedure $id $m $sub_group_open $id_list
		puts "$id\t$m\t$round_2_counter"
		close $sub_group_open
		set group_stat_array($m) "$group_members_counter"
		#### CREATE SUB-MATRIX FILE ####
		foreach sub_id1 $sub_group_list {
		    set sub_len [llength $sub_group_list]
		    if {$sub_len >= 2 && $sub_id1 != ""} {
			while {$sub_len >= 0} {
			    set sub_id2 [lindex $sub_group_list [expr $sub_len -1 ]]
			    if {$sub_id1 != $sub_id2 && $sub_id2 != ""} {
				# puts "$sub_id1\t$sub_id2"
				#####
				set query1 [info exists matrix_array($sub_id1,$sub_id2)]
				set query2 [info exists matrix_array($sub_id2,$sub_id1)]
				if {$query1 == 1} {
				    set idnt $matrix_array($sub_id1,$sub_id2)
				    puts $sub_matrx_open "$sub_id1\t$sub_id2\t$idnt"
				    unset matrix_array($sub_id1,$sub_id2)
				}
				if {$query2 == 1} {
				    set idnt $matrix_array($sub_id2,$sub_id1)
				    puts $sub_matrx_open "$sub_id2\t$sub_id1\t$idnt"
				    unset matrix_array($sub_id2,$sub_id1)
				}
				######
			    }
			    incr sub_len -1
			}
		    } 
		}
		close $sub_matrx_open
		#### END OF SUB-MATRIX FILE ####
		incr round_2_counter
	    }
	}
	puts "  end of round 2 "
	puts "=== DFS Done! ==="
	puts ""
	puts "Ready to generate final output files"
	puts ""
	after 1000
	Finish_Dummy_Groups $argv
    }

    proc DFS_Procedure {first_id m sub_group_open id_list} {

	global adj_array
	global used_nodes_list
	global group_members_counter
	global round_2_counter
	global retarded_counter
	global sub_group_list
	global group_array
	global list_of_group_array

	set l 1
	set group_counter 1
	foreach id $id_list {
	    set current_adj_list [split $adj_array($id)]
	    set perfect_match [lsearch -exact $current_adj_list $first_id]
	    if {$perfect_match >= 0} {
		set current_adj_line_length [llength [split $adj_array($id) " "]]
		while {$current_adj_line_length >= 0} {
		    set current_node_ID [lindex [split $adj_array($id) " "] [expr $current_adj_line_length -1 ]]
		    set already_done [lsearch -exact $used_nodes_list $current_node_ID]
		    if {$already_done < 0 && $current_node_ID != ""} {
			set group_array($current_node_ID) "$current_node_ID\t$m"
			set list_of_group_array [lappend list_of_group_array $group_array($current_node_ID)]
			puts $sub_group_open "$current_node_ID"
			set sub_group_list [lappend sub_group_list $current_node_ID]
			puts "$current_node_ID\t$m\t$round_2_counter"
			incr group_members_counter
			incr round_2_counter
			Already_Done $current_node_ID
			DFS_Procedure $current_node_ID $m $sub_group_open $id_list
			incr retarded_counter
		    }
		    incr current_adj_line_length -1
		}
	    }
	incr l
	}
	puts "Escape from DFS ... $retarded_counter  times"
    }

    proc Already_Done {current_node_ID} {

	global used_nodes_list

        set already_done [lsearch -exact $used_nodes_list $current_node_ID]
	if {$already_done < 0} {
	    set used_nodes_list [lappend used_nodes_list "$current_node_ID"]
	}
	if {$already_done >= 0} {
	    continue
	}
    }

    proc Finish_Dummy_Groups {argv} {

	global adj_array
	global group_array
	global group_stat_array
	global list_of_group_array

	set f_out8 [open [lindex $argv 1].group_degree_info "w"]

	foreach line $list_of_group_array {
	    set line [split $line]
	    set id [lindex $line 0]
	    set group [lindex $line 1]
	    set degree $group_stat_array($group)
	    set adjacent [expr [llength [split $adj_array($id)]] -1 ]
	    set mark [lindex $line 2]
	    puts $f_out8 "$id\t$adjacent\t$degree\t$group\t$mark"
	}
	close $f_out8

	puts ""
	puts "+---------------------------------+"
	puts "|     THE END OF GRAPH SEARCH     |"
	puts "+---------------------------------+"

    }

#### MAIN BODY #####

puts "$argc arguments entered"

if {$argc != 6} {

    puts "Program usage:"
    puts "tcl_blast_parser_123.tcl \[input\], \[output\], \[Exp(20)\], \[Idnt(40)\], \[Overlap(100)\], \[MATRIX/GRAPH option\]"

} else {

    puts $argv
    Blast_Parsing $argv

}

####  THE END  ####
