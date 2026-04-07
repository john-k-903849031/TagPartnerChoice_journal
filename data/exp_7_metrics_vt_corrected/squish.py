treatments_1 = [0, 1, 2]
treatments_2 = [0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625]
header_prefix = "seed,tag_mut,tag_perm,metric,"

def loop(seed, t_1, t_2, t_3, in_file_prefix, out_file):
    for metric in map(str, t_1): 
        for b in range(0,30):
            for tp in map(str, t_2):
                s = str(seed)
                seed += 1
                in_file_name = in_file_prefix + "_METRIC"+metric+"_TM"+t_3+"_TP"+tp+"_data_SEED" + s + ".data"
                line_prefix = s + "," + t_3 + "," + tp + "," + metric +","
                try :
                    in_file = open(in_file_name, 'r')
                    next(in_file)
                    for line in in_file :
                        out_file.write(line_prefix + line)
                    in_file.close()
                except :
                    print("could not read", in_file_name)



def munge_file(out_file_name, out_file_header, in_file_prefix) :
    print("writing", out_file_name + "....")
    out_file = open(out_file_name, 'w')
    out_file.write(out_file_header + "\n")
    s = 800001
    tms = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0]
    for tm in map(str, tms):    
        loop(s, treatments_1, treatments_2, tm, in_file_prefix, out_file)
        s += 810
    out_file.close()



header = header_prefix + "update,mean_intval,max_intval,min_intval,count,uninfected_host_count,Hist_-1,Hist_-0.9,Hist_-0.8,Hist_-0.7,Hist_-0.6,Hist_-0.5,Hist_-0.4,Hist_-0.3,Hist_-0.2,Hist_-0.1,Hist_0.0,Hist_0.1,Hist_0.2,Hist_0.3,Hist_0.4,Hist_0.5,Hist_0.6,Hist_0.7,Hist_0.8,Hist_0.9"
munge_file("host_counts.dat", header, "HostVals")

header = header_prefix + "update,host_mean_repro_count,sym_mean_repro_count,host_towards_partner_rate,host_from_partner_rate,sym_towards_partner_rate,sym_from_partner_rate"
munge_file("repro_hist.dat", header, "ReproHist")

header = header_prefix + "update,Variance_Hist_0,Variance_Hist_1,Variance_Hist_2,Variance_Hist_3,Variance_Hist_4,Variance_Hist_5,Variance_Hist_6,Variance_Hist_7,Variance_Hist_8,Variance_Hist_9,Variance_Hist_10,Variance_Hist_11,Variance_Hist_12,Variance_Hist_13,Variance_Hist_14,Variance_Hist_15,Variance_Hist_16,Variance_Hist_17,Variance_Hist_18,Variance_Hist_19,Variance_Hist_20,Mean_Hist_-1,Mean_Hist_-0.9,Mean_Hist_-0.8,Mean_Hist_-0.7,Mean_Hist_-0.6,Mean_Hist_-0.5,Mean_Hist_-0.4,Mean_Hist_-0.3,Mean_Hist_-0.2,Mean_Hist_-0.1,Mean_Hist_0.0,Mean_Hist_0.1,Mean_Hist_0.2,Mean_Hist_0.3,Mean_Hist_0.4,Mean_Hist_0.5,Mean_Hist_0.6,Mean_Hist_0.7,Mean_Hist_0.8"
munge_file("sym_diversity.dat", header, "SymDiversity")

header = header_prefix + "update,mean_intval,max_intval,min_intval,count,Hist_-1,Hist_-0.9,Hist_-0.8,Hist_-0.7,Hist_-0.6,Hist_-0.5,Hist_-0.4,Hist_-0.3,Hist_-0.2,Hist_-0.1,Hist_0.0,Hist_0.1,Hist_0.2,Hist_0.3,Hist_0.4,Hist_0.5,Hist_0.6,Hist_0.7,Hist_0.8,Hist_0.9"
munge_file("sym_counts.dat", header, "SymVals")

header = header_prefix + "update,mean_tag_distance,host_tag_richness,host_tag_shannon,symbiont_tag_richness,symbiont_tag_shannon,tag_0.1,tag_0.2,tag_0.3,tag_0.4,tag_0.5,tag_0.6,tag_0.7,tag_0.8,tag_0.9,tag_1.0"
munge_file("tag_dists.dat", header, "TagDist")

header = header_prefix + "update,horiz_attempt_-1_-0.8,horiz_attempt_-0.8_-0.6,horiz_attempt_-0.6_0.4,horiz_attempt_-0.4_0.2,horiz_attempt_-0.2_0,horiz_attempt_0_0.2,horiz_attempt_0.2_0.4,horiz_attempt_0.4_0.6,horiz_attempt_0.6_0.8,horiz_attempt_0.8_1,horiz_success_-1_-0.8,horiz_success_-0.8_-0.6,horiz_success_-0.6_0.4,horiz_success_-0.4_0.2,horiz_success_-0.2_0,horiz_success_0_0.2,horiz_success_0.2_0.4,horiz_success_0.4_0.6,horiz_success_0.6_0.8,horiz_success_0.8_1,vert_attempt_-1_-0.8,vert_attempt_-0.8_-0.6,vert_attempt_-0.6_0.4,vert_attempt_-0.4_0.2,vert_attempt_-0.2_0,vert_attempt_0_0.2,vert_attempt_0.2_0.4,vert_attempt_0.4_0.6,vert_attempt_0.6_0.8,vert_attempt_0.8_1,vert_success_-1_-0.8,vert_success_-0.8_-0.6,vert_success_-0.6_0.4,vert_success_-0.4_0.2,vert_success_-0.2_0,vert_success_0_0.2,vert_success_0.2_0.4,vert_success_0.4_0.6,vert_success_0.6_0.8,vert_success_0.8_1,horiz_tagfail_-1_-0.8,horiz_tagfail_-0.8_-0.6,horiz_tagfail_-0.6_0.4,horiz_tagfail_-0.4_0.2,horiz_tagfail_-0.2_0,horiz_tagfail_0_0.2,horiz_tagfail_0.2_0.4,horiz_tagfail_0.4_0.6,horiz_tagfail_0.6_0.8,horiz_tagfail_0.8_1,horiz_sizefail_-1_-0.8,horiz_sizefail_-0.8_-0.6,horiz_sizefail_-0.6_0.4,horiz_sizefail_-0.4_0.2,horiz_sizefail_-0.2_0,horiz_sizefail_0_0.2,horiz_sizefail_0.2_0.4,horiz_sizefail_0.4_0.6,horiz_sizefail_0.6_0.8,horiz_sizefail_0.8_1"
munge_file("trans_counts.dat", header, "TransmissionRates")

header = header_prefix + "host_int,sym_int,host_repro_count,host_towards_partner_count,host_from_partner_count,sym_repro_count,sym_towards_partner_count,sym_from_partner_count,host_tag,sym_tag,tag_distance"
munge_file("org_dump.dat", header, "OrgDump")
