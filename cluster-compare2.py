#!/usr/bin/env python
# vim:foldmethod=marker
# combine and analyse HIV-Trace"s output

### IMPORTING MODULES {{{
import sys
import argparse
import numpy as np
import pandas as pd
import json as js
import statistics as st
from igraph import *
from subprocess import Popen
from collections import Counter
from pathlib import Path
from datetime import datetime, timedelta
from scipy import stats
from seaborn import color_palette
from random import uniform
#}}}

### PATHS & ARGUMENTS {{{
settings_file = "hivtrace_settings.csv"         # contains settings used with which HIV-Trace was run
cl_stats_file = "cc_cl_stats.csv"               # contains cluster-statistics of all timepoints
cl_comp_file = "cc_cl_comp.csv"                 # contains all the comparisons and scores between clusters in subsequent timepoints
cl_traj_file = "cc_cl_traj.csv"                 # contains the trajectories of each unique cluster
cl_size_file = "cc_cl_size.csv"                 # contains size-trajectories of each unique cluster
cl_class_file = "cc_cl_class.csv"               # contains cluster classification along chosen variable(s)
nd_stats_file = "cc_nd_stats.csv"               # contains node-statistics of all timepoints
nd_stats_bl_file = "cc_bl_nd_stats.csv"         # contains node-statistics at baseline
nd_attrib_file = "cc_attributes.csv"            # contains combined statistics of each node at each timepoint 
nd_attrib_bl_file = "cc_bl_attributes.csv"      # contains combined statistics of each node at baseline
cl_v_file = "cl_v.csv"                          # contains predictions and outcomes of a cluster that's plotted using --pred_plot
tables_dir = "/home/cluster/mlabar/data/2103stata/"      # path to stata-tables containing patient-data


parser = argparse.ArgumentParser(description="Follow clusters through multiple HIV-Trace output files.")

parser.add_argument("-f", "--fasta", help = "Input FASTA file", dest = "fasta", type = str, required = True)
parser.add_argument("-j", "--json", help = "Input JSON file", dest = "json", type = str, required = True)
parser.add_argument("-d", "--date", help = "Input date", dest = "dates", type = str, required = True)
parser.add_argument("-p", "--patient_id", help = "Patient-ID of the new patient", dest = "pid", type = str, required = True)
parser.add_argument("--skip_baseline", help = "Do not re-calculate node- and cluster statistics at baseline", dest = "skip_baseline", action="store_true")
parser.add_argument("--outcome_interval", help = "The time interval over which an outcome is caluclated in years", dest = "outcome_interval", type = int, default = 3)
parser.add_argument("--classify_clusters", help = "Classify clusters according to one or many variables", dest = "classify", nargs = "+")
parser.add_argument("-v", help = "Set verbosity", dest = "verbosity", type = int, choices = [0, 1, 2], default = 2)

args = parser.parse_args()

fasta = args.fasta
json = args.json
date = args.dates
new_patient = int(args.pid)
#}}}

### FUNCTIONS {{{
def get_settings(json):
    with open(json) as c:
        jsl = js.load(c)
    
    thresh = [jsl["trace_results"]["Settings"]["threshold"]]
    edge_filter = [jsl["trace_results"]["Settings"]["edge-filtering"]]
    cont = [jsl["trace_results"]["Settings"]["contaminants"]]
    singletons = [jsl["trace_results"]["Settings"]["singletons"]]
    compact_json = [jsl["trace_results"]["Settings"]["compact_json"]]
    created = [jsl["trace_results"]["Settings"]["created"]]
    set_tbl = pd.DataFrame(data =   {"threshold": thresh, 
                                    "edge_filtering": edge_filter, 
                                    "contaminants": cont, 
                                    "singletons": singletons, 
                                    "compact_json": compact_json, 
                                    "creation date": created})
    
    return set_tbl

def get_clusters(json):
    with open(json) as c:
        cl_sizes = js.load(c)
        
    cl_sizes = cl_sizes["trace_results"]["Cluster sizes"]
    
    cl_sizes.sort(reverse = True) # by default, the list of cluster sizes is unsorted! (why??)
    
    cl_tbl = pd.DataFrame(data = {"cluster": range(1, len(cl_sizes) + 1), "size": cl_sizes})

    return cl_tbl

def get_nodes(json):
    with open(json) as c:
        cl_nodes = js.load(c)
    
    ids = cl_nodes["trace_results"]["Nodes"]["id"]
    clusters = cl_nodes["trace_results"]["Nodes"]["cluster"]
    
    try:
        nd_tbl = pd.DataFrame(data = {"id": ids, "cluster": clusters})
    except ValueError:
        nd_val_lst = clusters["values"]
        nd_val_lst_2 = [x + 1 for x in nd_val_lst]
        nd_dict = {"id": ids, "cluster": nd_val_lst_2}
        nd_tbl = pd.DataFrame(data = nd_dict)

    return nd_tbl

def get_edges(json):
    with open(json) as c:
        cl_edges = js.load(c)
    
    ed_seq = cl_edges["trace_results"]["Edges"]["sequences"]
    ed_len = cl_edges["trace_results"]["Edges"]["length"]
    ed_src = cl_edges["trace_results"]["Edges"]["source"]
    ed_tgt = cl_edges["trace_results"]["Edges"]["target"]

    ed_seq1 = []
    ed_seq2 = []
    for pair in ed_seq:
        ed_seq1.append(pair[0])
        ed_seq2.append(pair[1])

    try:
        ed_tbl = pd.DataFrame(data = {"seq1": ed_seq1, "seq2": ed_seq2, "length": ed_len})
    except ValueError:
        ed_l_keys = ed_len["keys"]
        ed_l_values = ed_len["values"]

        ed_len2 = []
        for val in ed_l_values:
            ed_len2.append(ed_l_keys[str(val)])

        ed_tbl = pd.DataFrame(data = {"seq1": ed_seq1, "seq2": ed_seq2, "length": ed_len2})
    
    return ed_tbl

def cluster_graph(cluster_nr, nd_tbl, ed_tbl):
    nw_nodes = nd_tbl[nd_tbl["cluster"] == int(cluster_nr)]
    nw_nodes_lst = nw_nodes["id"].tolist()

    # collect edges
    ed1_mask = [(x in nw_nodes_lst) for x in ed_tbl.seq1.tolist()]
    ed2_mask = [(x in nw_nodes_lst) for x in ed_tbl.seq2.tolist()]
    ed_mask = [(ed1_mask[i] == True | ed2_mask[i] == True) for i in range(len(ed1_mask))]

    ed_tbl2 = ed_tbl.loc[ed_mask, ]

    ed_tbl2 = ed_tbl2.drop_duplicates()
    
    # the old way: writing edge-list to file and reading it with Graph.Read_Ncol
    ed_tbl2.to_csv("readjson_edges.txt", index = False, header = False, sep = "\t") 
    G = Graph.Read_Ncol("readjson_edges.txt", directed = False)
    
    #ed_tbl2 = ed_tbl2.rename(columns = {"length": "weight"})
    #G = Graph.DataFrame(ed_tbl2, directed = False, vertices = nw_nodes)

    return G

def cluster_stats(graph):
    if graph == 0:
        node_n = 0
        deg_range = 0
        deg_iqr = 0
        deg_mean = 0
        deg_median = 0
        dist_range = 0
        dist_iqr = 0
        dist_mean = 0
        dist_median = 0
        betw_median = 0
        close_median = 0
        dens = 0
        diam = 0
        trans = 0
    
    else:
        degrees = []
        for i in range(len(graph.vs)):
            degrees.append(graph.degree(i))
        
        node_n = len(graph.vs)
        distances = graph.es["weight"]

        deg_mean = st.mean(degrees)
        deg_median = st.median(degrees)
        deg_range = max(degrees) - min(degrees)
        deg_iqr = stats.iqr(degrees)

        dist_mean = st.mean(distances)
        dist_median = st.median(distances)
        dist_median = distances[0]
        dist_range = max(distances) - min(distances)
        dist_iqr = stats.iqr(distances)

        betw_median = st.median(graph.betweenness())
        close_median = st.median(graph.closeness())
        dens = graph.density()
        diam = graph.diameter(weights = "weight")
        trans = graph.transitivity_undirected()

    cl_stats = {"cl_size": node_n, 
                "cl_degree_range": deg_range, 
                "cl_degree_iqr": deg_iqr, 
                "cl_degree_mean": deg_mean, 
                "cl_degree_median": deg_median, 
                "cl_distance_range": dist_range, 
                "cl_distance_iqr": dist_iqr, 
                "cl_distance_mean": dist_mean, 
                "cl_distance_median": dist_median,
                "cl_betweenness_median": betw_median,
                "cl_closeness_median": close_median,
                "cl_density": dens,
                "cl_diameter": diam,
                "cl_transitivity": trans}

    return cl_stats

def calculate_present(new_patient_cluster, nd_tbl, ed_tbl):
    # calculate cluster-statistics of nodes at present
    if new_patient_cluster != 0:
        G = cluster_graph(new_patient_cluster, nd_tbl, ed_tbl)
        name_list = G.vs["name"]
        pid_list = [x.split("|")[1] for x in name_list]
        timestamp_list = [x.split("|")[2] for x in name_list]
        v_G = pid_list.index(str(new_patient))
        node_n_G = len(G.vs) # present size of the cluster

        # calculating past size for cluster growth
        past_timestamp = datetime.strftime((datetime.strptime(new_patient_timestamp, "%Y-%m-%d %H:%M:%S") - timedelta(days = 365 * args.outcome_interval)), "%Y-%m-%dT%H:%M:%S")
        vert_mask = [x > past_timestamp for x in timestamp_list]
        node_n_G_past = vert_mask.count(False)

        cl_gr_past = node_n_G - node_n_G_past

        degrees = G.degree()
        distances = G.es["weight"]
        betws = G.betweenness()
        close = G.closeness()
        diam = G.diameter(weights = "weight")

        deg_mean = st.mean(degrees)
        deg_median = st.median(degrees)
        deg_range = max(degrees) - min(degrees)
        deg_iqr = stats.iqr(degrees)

        nd_deg = degrees[v_G]
        nd_close = close[v_G]
        dist_mean = st.mean(distances)
        dist_median = distances[0]
        betw_median = st.median(betws)
        close_median = st.median(close)
        dist_range = max(distances) - min(distances)
        dist_iqr = stats.iqr(distances)
        dens = G.density()
        trans = G.transitivity_undirected()
        
        cl_stats = {"id": new_patient_header, 
                    "patient_id": new_patient,
                    "sequence_timestamp": new_patient_timestamp,
                    "cl_size": node_n_G, 
                    "cl_past_gr": cl_gr_past, 
                    "cl_degree_range": deg_range, 
                    "cl_degree_iqr": deg_iqr, 
                    "cl_degree_mean": deg_mean, 
                    "cl_degree_median": deg_median, 
                    "cl_distance_range": dist_range, 
                    "cl_distance_iqr": dist_iqr, 
                    "cl_distance_mean": dist_mean, 
                    "cl_distance_median": dist_median,
                    "cl_betweenness_median": betw_median,
                    "cl_closeness_median": close_median,
                    "cl_density": dens,
                    "cl_diameter": diam,
                    "cl_transitivity": trans, 
                    "nd_degree": nd_deg, 
                    "nd_betweenness": betws[v_G], 
                    "nd_closeness": nd_close}
    else:
        cl_stats = {"id": new_patient_header, 
                    "patient_id": new_patient,
                    "sequence_timestamp": new_patient_timestamp,
                    "cl_size": 1, 
                    "cl_past_gr": 1, 
                    "cl_gr": np.nan, 
                    "cl_degree_range": 0, 
                    "cl_degree_iqr": 0, 
                    "cl_degree_mean": 0, 
                    "cl_degree_median": 0, 
                    "cl_distance_range": 0, 
                    "cl_distance_iqr": 0, 
                    "cl_distance_mean": 0, 
                    "cl_distance_median": 0,
                    "cl_betweenness_median": 0,
                    "cl_closeness_median": 0,
                    "cl_density": 0,
                    "cl_diameter": 0,
                    "cl_transitivity": 0, 
                    "nd_degree": 0, 
                    "nd_gr": np.nan, 
                    "nd_betweenness": 0, 
                    "nd_closeness": 0}

    cl_stats_df = pd.DataFrame(cl_stats, index = [0])

    return cl_stats_df

def calculate_baseline(header_list, cl_list, nd_tbl, ed_tbl, ignore_list = []):
    # calculate cluster-statistics of nodes at baseline
    # process clustered sequences
    both_bl_list = []
    stats_bl_list = []
    for cl in cl_list:
        G = cluster_graph(cl, nd_tbl, ed_tbl)
        name_list = G.vs["name"]
        pid_list = [x.split("|")[1] for x in name_list]
        timestamp_list = [x.split("|")[2] for x in name_list]

        for v in range(len(G.vs)):
            #print(">>>Cluster: " + str(cl) + " (" + str(v) + " / " + str(len(G.vs)) + ")", end = "\r")

            v_name = name_list[v]
            v_pid = pid_list[v]
            v_timestamp = timestamp_list[v]

            if v_name in ignore_list:
                continue

            past_timestamp = datetime.strftime((datetime.strptime(v_timestamp, "%Y-%m-%dT%H:%M:%S") - timedelta(days = 365 * args.outcome_interval)), "%Y-%m-%dT%H:%M:%S")
            future_timestamp = datetime.strftime((datetime.strptime(v_timestamp, "%Y-%m-%dT%H:%M:%S") + timedelta(days = 365 * args.outcome_interval)), "%Y-%m-%dT%H:%M:%S")

            # make a subgraph containing all the nodes that are older (or equally old) than node v
            vert_mask = [x > timestamp_list[v] for x in timestamp_list] # boolean mask that is True for all vertex ids that have a newer timestamp than the vertex with id = v
            
            del_vert_ids = []
            for i, bool_ in enumerate(vert_mask):
                if bool_:
                    del_vert_ids.append(i)

            Gsub = G.copy()
            Gsub.delete_vertices(del_vert_ids) # new graph with all newer sequences than vertex v removed
            v_Gsub = Gsub.vs.find(name = v_name).index # new vertex-ID within Gsub
            
            comps = list(Gsub.clusters()) # list of connected components
            for i, comp in enumerate(comps):
                if v_Gsub in comp:
                    Gsub_conn = Gsub.clusters().subgraph(i).copy() # new graph of only the connected component that contains vertex with ID v_Gsub
                    v_Gsub_conn = Gsub_conn.vs.find(name = v_name).index # another new vertex-ID within Gsub_conn

                    # make a subset of the above graph containing only nodes equal or older than the past timestamp (necessary for calculation of cluster growth)
                    Gsub_conn_timestamps = [x.split("|")[2] for x in list(Gsub_conn.vs["name"])]
                    vert_mask_conn = [x > past_timestamp for x in Gsub_conn_timestamps]

                    del_vert_ids_conn = []
                    for j, bool_ in enumerate(vert_mask_conn):
                        if bool_:
                            del_vert_ids_conn.append(j)

                    Gsub_conn_sub = Gsub_conn.copy()
                    Gsub_conn_sub.delete_vertices(del_vert_ids_conn)

                    node_n_conn_sub = len(Gsub_conn_sub.vs)

                    break

            node_n_conn = len(Gsub_conn.vs)
            node_n_conn_sub = len(Gsub_conn_sub.vs)
            nd_deg = Gsub_conn.degree(v_Gsub_conn)

            # make another subgraph to calculate future cluster and node growth
            if future_timestamp <= date:
                vert_mask_3 = [x > future_timestamp for x in timestamp_list]

                del_vert_ids_3 = []
                for i, bool_ in enumerate(vert_mask_3):
                    if bool_:
                        del_vert_ids_3.append(i)

                Gsub_3 = G.copy()
                Gsub_3.delete_vertices(del_vert_ids_3)
                v_Gsub_3 = Gsub_3.vs.find(name = v_name).index # new vertex-ID within Gsub_3
            
                comps_3 = list(Gsub_3.clusters()) # list of connected components
                for i, comp in enumerate(comps_3):
                    if v_Gsub_3 in comp:
                        Gsub_3_conn = Gsub_3.clusters().subgraph(i).copy() # new graph of only the connected component that contains vertex with ID v_Gsub
                        v_Gsub_3_conn = Gsub_3_conn.vs.find(name = v_name).index # another new vertex-ID within Gsub_conn
                        
                        # make a subset of the above graph containing only nodes equal or older than the node v (necessary for calculation of cluster growth)
                        Gsub_3_conn_timestamps = [x.split("|")[2] for x in list(Gsub_3_conn.vs["name"])]
                        vert_mask_3_conn = [x > v_timestamp for x in Gsub_3_conn_timestamps]

                        del_vert_ids_3_conn = []
                        for j, bool_ in enumerate(vert_mask_3_conn):
                            if bool_:
                                del_vert_ids_3_conn.append(j)

                        Gsub_3_conn_sub = Gsub_3_conn.copy()
                        Gsub_3_conn_sub.delete_vertices(del_vert_ids_3_conn)

                        node_n_new_conn_sub = len(Gsub_3_conn_sub.vs)

                        break

                node_n_new_conn = len(Gsub_3_conn.vs)
                node_n_new_conn_sub = len(Gsub_3_conn_sub.vs)
                nd_deg_new = Gsub_3_conn.degree(v_Gsub_3_conn)

            else:
                node_n_new_conn = np.nan
                node_n_new_conn_sub = np.nan
                nd_deg_new = np.nan


            # calculating statistics
            cl_gr_past = node_n_conn - node_n_conn_sub
            cl_gr = node_n_new_conn - node_n_new_conn_sub
            nd_gr = nd_deg_new - nd_deg
            if nd_gr > 0:
                nd_gr_bool = True
            elif nd_gr == 0:
                nd_gr_bool = False
            else:
                nd_gr_bool = np.nan

            degrees = Gsub_conn.degree()
            distances = Gsub_conn.es["weight"]
            betws = Gsub_conn.betweenness()
            close = Gsub_conn.closeness()
            diam = Gsub_conn.diameter(weights = "weight")

            deg_mean = st.mean(degrees)
            deg_median = st.median(degrees)
            deg_range = max(degrees) - min(degrees)
            deg_iqr = stats.iqr(degrees)

            if len(distances) > 0:
                dist_mean = st.mean(distances)
                dist_median = distances[0]
                betw_median = st.median(betws)
                close_median = st.median(close)
                dist_range = max(distances) - min(distances)
                dist_iqr = stats.iqr(distances)
                dens = Gsub_conn.density()
                nd_close = close[v_Gsub_conn]
            else:
                dist_mean = 0
                dist_median = 0
                betw_median = 0
                close_median = 0
                dist_range = 0
                dist_iqr = 0
                dens = 0
                nd_close = 0

            if len(distances) > 1:
                trans = Gsub_conn.transitivity_undirected()
            else:
                trans = 0
            
            cl_bl_stats = {"id": v_name, 
                        "patient_id": v_pid,
                        "sequence_timestamp": v_timestamp,
                        "cl_size": node_n_conn, 
                        "cl_past_gr": cl_gr_past, 
                        "cl_gr": cl_gr, 
                        "cl_degree_range": deg_range, 
                        "cl_degree_iqr": deg_iqr, 
                        "cl_degree_mean": deg_mean, 
                        "cl_degree_median": deg_median, 
                        "cl_distance_range": dist_range, 
                        "cl_distance_iqr": dist_iqr, 
                        "cl_distance_mean": dist_mean, 
                        "cl_distance_median": dist_median,
                        "cl_betweenness_median": betw_median,
                        "cl_closeness_median": close_median,
                        "cl_density": dens,
                        "cl_diameter": diam,
                        "cl_transitivity": trans, 
                        "nd_degree": nd_deg, 
                        "nd_gr": nd_gr, 
                        "nd_gr_bool": nd_gr_bool,
                        "nd_betweenness": betws[v_Gsub_conn], 
                        "nd_closeness": nd_close}

            stats_bl_df = pd.DataFrame(cl_bl_stats, index = [0])
            stats_bl_list.append(stats_bl_df)

    if len(stats_bl_list) > 1:
        stats_bl = pd.concat(stats_bl_list)
        both_bl_list.append(stats_bl)
    elif len(stats_bl_list) == 1:
        stats_bl = stats_bl_list[0]
        both_bl_list.append(stats_bl)

    # add singletons
    try:
        clustered_ids = stats_bl["id"].tolist()
    except:
        clustered_ids = []

    singleton_stats_bl_list = []
    for id_ in header_list:
        if id_ in clustered_ids or id_ in ignore_list:
            continue

        v_pid, v_timestamp = id_.split("|")[1:]

        future_timestamp = datetime.strftime((datetime.strptime(v_timestamp, "%Y-%m-%dT%H:%M:%S") + timedelta(days = 365 * args.outcome_interval)), "%Y-%m-%dT%H:%M:%S")
        if future_timestamp > date:
            nd_gr = np.nan
            nd_gr_bool = np.nan
        else:
            nd_gr = 0
            nd_gr_bool = False

        cl_bl_stats = {"id": id_, 
                "patient_id": v_pid,
                "sequence_timestamp": v_timestamp,
                "cl_size": 1, 
                "cl_past_gr": 1, 
                "cl_gr": 0, 
                "cl_degree_range": 0, 
                "cl_degree_iqr": 0, 
                "cl_degree_mean": 0, 
                "cl_degree_median": 0, 
                "cl_distance_range": 0, 
                "cl_distance_iqr": 0, 
                "cl_distance_mean": 0, 
                "cl_distance_median": 0,
                "cl_betweenness_median": 0,
                "cl_closeness_median": 0,
                "cl_density": 0,
                "cl_diameter": 0,
                "cl_transitivity": 0, 
                "nd_degree": 0, 
                "nd_gr": nd_gr, 
                "nd_gr_bool": nd_gr_bool,
                "nd_betweenness": 0, 
                "nd_closeness": 0}
    
        cl_bl_stats_df = pd.DataFrame(cl_bl_stats, index = [0])
        singleton_stats_bl_list.append(cl_bl_stats_df)

    if len(singleton_stats_bl_list) > 1:
        singleton_stats_bl = pd.concat(singleton_stats_bl_list)
        both_bl_list.append(singleton_stats_bl)
    elif len(singleton_stats_bl_list) == 1:
        singleton_stats_bl = singleton_stats_bl_list[0]
        both_bl_list.append(singleton_stats_bl)

    if len(both_bl_list) > 1:
        stats_bl = pd.concat(both_bl_list)
        stats_bl = stats_bl.astype({"patient_id": "int32"})
    elif len(both_bl_list) == 1:
        stats_bl = both_bl_list[0]
        stats_bl = stats_bl.astype({"patient_id": "int32"})
    else:
        stats_bl = pd.DataFrame()

    return stats_bl

def filter_table_present(patient_id, date, tbl, datevar, tbl_vars, date_format = "%Y-%m-%d"):
    tbl = tbl.replace("", np.nan)                                                           # replace empty strings in the date-column with NAs
    tbl[datevar] = pd.to_datetime(tbl[datevar], format = date_format, errors = "coerce")    # cast to datetime and return NaT for dates that produce errors
    tbl = tbl.dropna(subset = [datevar])                                                    # remove rows without valid date

    subtbl = tbl[tbl["patient_id"] == patient_id]
    subtbl = subtbl.sort_values(by = datevar, ascending = False)
    subtbl = subtbl[subtbl[datevar] <= date] # remove rows older than date

    if subtbl.shape[0] > 0:
        tbl_filtered = subtbl.iloc[[0]] # newest row according to datevar
    else:
        tmp_dict_1 = {"patient_id": patient_id}
        tmp_dict_2 = dict.fromkeys(tbl_vars, np.NaN)
        tmp_dict = {**tmp_dict_1, **tmp_dict_2} # combine the above dicts into one

        tbl_filtered = pd.DataFrame(tmp_dict, index = [0])

    tbl_filtered = tbl_filtered.astype({"patient_id": "int32"})
    
    return tbl_filtered

def filter_table_baseline(tbl, datevar, pid_to_ts_dict, date_format = "%Y-%m-%d"): # same as filter_table, but return a table with the very first observation of each patient
    tbl = tbl.replace("", np.nan)                                                           # replace empty strings in the date-column with NAs
    tbl[datevar] = pd.to_datetime(tbl[datevar], format = date_format, errors = "coerce")    # cast to datetime and return NaT for dates that produce errors
    tbl = tbl.dropna(subset = [datevar])                                                    # remove rows without valid date

    u_pids_tbl = list(set(tbl["patient_id"].tolist()))
    u_pids = []
    for u_pid in u_pids_tbl:
        if u_pid in cl_stats_bl["patient_id"].tolist(): # keep only the patient-ids that have associated degree-data (skips some time-consuming filtering in the next step)
            u_pids.append(u_pid)

    row_list = []
    for u_pid in u_pids:
        ts = pid_to_ts_dict[u_pid]

        subtbl = tbl[tbl["patient_id"] == u_pid]                        # filter so only one patient's entries remain
        subtbl = subtbl.sort_values(by = datevar, ascending = True)    # sort rows by datevar with the most recent row on the bottom
        
        subtbl_dates = subtbl[datevar].tolist()
        subtbl_date_deltas = [abs((ts - x).days) for x in subtbl_dates] # list of timespans between datevar and sequence timestamp
        index = subtbl_date_deltas.index(min(subtbl_date_deltas))   # index of the entry that is closest to the sequence timestamp

        subtbl = subtbl.iloc[[index]]

        row_list.append(subtbl)

    tbl_filtered = pd.concat(row_list)
    tbl_filtered = tbl_filtered.astype({"patient_id": "int32"})
    
    return tbl_filtered

def classify_clusters(var_list):
    all_steps = []
    for step in range(steps):
        cl_nrs = cl_list[step]
        step_rows = []
        for cl_nr in cl_nrs:
            cl_uid = traj_dict_inv(step, cl_nr, traj_dict)

            cl_df = degs_attrib_df[degs_attrib_df["nd_cluster_" + str(step)] == cl_uid] # all nodes that are in the cluster at timepoint %step%
            cl_stats_df = cl_stats_tbl[(cl_stats_tbl["cluster_uid"] == cl_uid) & (cl_stats_tbl["step"] == step)] # cluster-stats of the cluster at timepoint %step%
            cl_size = cl_stats_df["cl_size"].tolist()[0]

            row = pd.DataFrame({"cluster": [cl_uid], 
                                "step": [step], 
                                "cl_size": [cl_size]})

            for var in var_list:
                if var[-1] == "_": # account for variables that are recorded once for each timestep
                    var = var + str(step)

                vals = cl_df[var].tolist()

                val_counts = Counter(vals)
                val_most_common = val_counts.most_common()

                # recording the most common value
                mc_val, mc_count = val_most_common[0]
                row["most_common_" + var] = [mc_val]
                frac = mc_count / len(vals)
                row["most_common_" + var + "_fraction"] = [frac]

                # recording the second-most common value
                if len(val_most_common) > 1:
                    mc_val, mc_count = val_most_common[1]
                    frac = mc_count / len(vals)
                else:
                    mc_val, mc_count = (np.nan, np.nan)
                    frac = np.nan

                row["second_most_common_" + var] = [mc_val]
                row["second_most_common_" + var + "_fraction"] = [frac]
                    
                # recording the third-most common value
                if len(val_most_common) > 2:
                    mc_val, mc_count = val_most_common[2]
                    frac = mc_count / len(vals)
                else:
                    mc_val, mc_count = (np.nan, np.nan)
                    frac = np.nan

                row["third_most_common_" + var] = [mc_val]
                row["third_most_common_" + var + "_fraction"] = [frac]

            step_rows.append(row)

        step_rows = pd.concat(step_rows)
        all_steps.append(step_rows)

    all_steps = pd.concat(all_steps)

    return all_steps
#}}}

### SETUP {{{
# make sure settings-file exists so the used threshold can be read
print("Generating and reading " + str(settings_file))
set_tbl = get_settings(json)
set_tbl.to_csv(settings_file, index = False)

hivtrace_settings = pd.read_csv(settings_file)
hivtrace_threshold = hivtrace_settings.iloc[0, 0]

# Collect headers from FASTA-files
print("Collecting headers from FASTA files")
with open(fasta) as f:
    lines = f.readlines()

header_list = []
pid_to_timestamp = {}
for i, line in enumerate(lines):
    if line[0] == ">":
        h = lines[i].strip().strip(">")
        header_list.append(h)

        pid = int(h.split("|")[1])
        timestamp = h.split("|")[2].split("T")[0]
        pid_to_timestamp[pid] = datetime.strptime(timestamp, "%Y-%m-%d")

        if int(h.split("|")[1]) == new_patient:
            new_patient_header = h
            new_patient_timestamp = datetime.strftime(datetime.strptime(timestamp, "%Y-%m-%d"), "%Y-%m-%d %H:%M:%S")

# getting node and edge lists for each JSON-file
print("Collecting node and edge lists")
nd_tbl = get_nodes(json)
ed_tbl = get_edges(json)

nd_tbl_ids = nd_tbl["id"].tolist()
nd_tbl_pids = [x.split("|")[1] for x in nd_tbl_ids]

try:
    nd_tbl_mask = nd_tbl_pids.index(str(new_patient))
    new_patient_cluster = nd_tbl["cluster"].tolist()[nd_tbl_mask]
except ValueError:
    new_patient_cluster = 0

# get cluster lists
print("Reading cluster lists")
cl_tbl = get_clusters(json)
cl_list = cl_tbl["cluster"].tolist()
#}}}

### ANALYSIS {{{
# Calculating cluster statistics at present (only the cluster that new_patient is in)
print("Calculating cluster statistics and exporting to " + str(cl_stats_file))
cl_stats = calculate_present(new_patient_cluster, nd_tbl, ed_tbl)
cl_stats.to_csv(cl_stats_file, index = False)

print("Calculating cluster statistics at baseline and exporting to " + str(nd_stats_bl_file))

if not args.skip_baseline:
    if Path(nd_stats_bl_file).is_file():
        cl_stats_bl_old = pd.read_csv(nd_stats_bl_file)
        outcome_timestamp = datetime.strftime((datetime.strptime(date, "%Y-%m-%d") - timedelta(days = 365 * args.outcome_interval)), "%Y-%m-%dT%H:%M:%S")

        ign_list1 = cl_stats_bl_old[cl_stats_bl_old["sequence_timestamp"] > outcome_timestamp]["id"].tolist() # sequence-IDs of patients that have been sampled <3 years ago and therefore have no outcome
        ign_list2 = cl_stats_bl_old[~cl_stats_bl_old["nd_gr"].isna()]["id"].tolist() # sequence-IDs of patients that DO have an outcome and therefore do not have to be recalculated
        ign_list = ign_list1 + ign_list2 # the sequence-IDs that have not changed since the last iteration and are therefore ignored

        ign_list_complement = [] # the sequence-IDs that have changed and are therefore recalculated with calculate_baseline()
        for id_ in cl_stats_bl_old["id"].tolist():
                if id_ not in ign_list:
                        ign_list_complement.append(id_)
    
        cl_stats_bl_new = calculate_baseline(header_list, cl_list, nd_tbl, ed_tbl, ign_list) # the new part of the cluster stats
        cl_stats_bl_old_mask = ~cl_stats_bl_old["id"].isin(ign_list_complement)
        cl_stats_bl_old_unchanged = cl_stats_bl_old[cl_stats_bl_old_mask] # the unchanged part of the cluster stats that is simply copied from the previous version

        cl_stats_bl = pd.concat([cl_stats_bl_old_unchanged, cl_stats_bl_new])
        cl_stats_bl.to_csv(nd_stats_bl_file, index = False)

    else:
        cl_stats_bl = calculate_baseline(header_list, cl_list, nd_tbl, ed_tbl)
        cl_stats_bl.to_csv(nd_stats_bl_file, index = False)

else:
    cl_stats_bl = pd.read_csv(nd_stats_bl_file)
#}}}

### PATIENT DATA IMPORT {{{
# Defining the variables to be imported for each table
pat_tbl_vars = ["id", "born", "sex", "center1", "last_center", "ethnicity", "hiv_posdate", "x_pref", "risk", "education", "profession", "virus_type", "infect_place"]
lab_tbl_vars = ["id", "labdate", "sample_day", "leu", "lym", "cd3", "cd4", "cd8", "rna", "syph_date", "syph_q", "hbvdna_qual", "antihcv"]
fup_tbl_vars = ["id", "fupdate", "p_stable_sex", "p_stable_cuse", "p_occas", "p_occas_sex", "p_occas_cuse", "p_occas_osex", "p_occas_number", "alcohol", "alc_freq", "hero_iv", "hero_iv_f", "coca_iv", "coca_iv_f", "other_iv", "other_iv_f", "hero_ni", "hero_ni_f", "coca_ni", "coca_ni_f", "cana_ni", "cana_ni_f", "other_ni", "other_ni_f"]
drug_tbl_vars = ["id", "drug", "starts", "stops", "stop_why"]
adhe_tbl_vars = ["id", "ad_date", "missed", "in_row"]

# Reading patient data from stata-tables
print("Reading patient data")
pat_tbl = pd.read_stata(tables_dir + "pat.dta", preserve_dtypes = False, columns = pat_tbl_vars)
lab_tbl = pd.read_stata(tables_dir + "lab.dta", preserve_dtypes = False, columns = lab_tbl_vars)
fup_tbl = pd.read_stata(tables_dir + "fup.dta", preserve_dtypes = False, columns = fup_tbl_vars)
drug_tbl = pd.read_stata(tables_dir + "drug.dta", preserve_dtypes = False, columns = drug_tbl_vars)
adhe_tbl = pd.read_stata(tables_dir + "adhe.dta", preserve_dtypes = False, columns = adhe_tbl_vars)

for tbl in [pat_tbl, lab_tbl, fup_tbl, drug_tbl, adhe_tbl]:
    tbl.rename(columns = {"id": "patient_id"}, inplace = True)

# select row of the new patient
pat_tbl = pat_tbl.astype({"patient_id": "int32"})

# Filtering tables to include the most recent row in each step
print("Filtering patient data")
filter_date = datetime.strftime((datetime.strptime(date, "%Y-%m-%d") + timedelta(days = 30)), "%Y-%m-%d")

lab_filtered = filter_table_present(new_patient, filter_date, lab_tbl, "sample_day", lab_tbl_vars[1: ], date_format = "%d/%m/%Y")
fup_filtered = filter_table_present(new_patient, filter_date, fup_tbl, "fupdate", fup_tbl_vars[1: ])
adhe_filtered = filter_table_present(new_patient, filter_date, adhe_tbl, "ad_date", adhe_tbl_vars[1: ])

# prepending the table of origin for each variable
pat_tbl = pat_tbl.rename(lambda a: "pat_" + a, axis = 1)
pat_tbl = pat_tbl.rename(columns = {"pat_patient_id": "patient_id"}) # rename patient_id-column back to the original name for easier merging later

lab_filtered = lab_filtered.rename(lambda a: "lab_" + a, axis = 1)
lab_filtered = lab_filtered.rename(columns = {"lab_patient_id": "patient_id"})

fup_filtered = fup_filtered.rename(lambda a: "fup_" + a, axis = 1)
fup_filtered = fup_filtered.rename(columns = {"fup_patient_id": "patient_id"})

adhe_filtered = adhe_filtered.rename(lambda a: "adhe_" + a, axis = 1)
adhe_filtered = adhe_filtered.rename(columns = {"adhe_patient_id": "patient_id"})

# Same for baseline values of each patient
# Filtering tables to include only the baseline-value for each patient
print("Filtering patient data at baseline")
lab_filtered_bl = filter_table_baseline(lab_tbl, "sample_day", pid_to_timestamp, date_format = "%d/%m/%Y")
fup_filtered_bl = filter_table_baseline(fup_tbl, "fupdate", pid_to_timestamp)
adhe_filtered_bl = filter_table_baseline(adhe_tbl, "ad_date", pid_to_timestamp)

# prepending the table of origin for each variable
lab_filtered_bl = lab_filtered_bl.rename(lambda a: "lab_" + a, axis = 1)
lab_filtered_bl = lab_filtered_bl.rename(columns = {"lab_patient_id": "patient_id"})

fup_filtered_bl = fup_filtered_bl.rename(lambda a: "fup_" + a, axis = 1)
fup_filtered_bl = fup_filtered_bl.rename(columns = {"fup_patient_id": "patient_id"})

adhe_filtered_bl = adhe_filtered_bl.rename(lambda a: "adhe_" + a, axis = 1)
adhe_filtered_bl = adhe_filtered_bl.rename(columns = {"adhe_patient_id": "patient_id"})

# Merging tables
print("Merging attribute tables")
pat_tbl_np = pat_tbl[pat_tbl["patient_id"] == new_patient]

attrib_df = pd.merge(pat_tbl_np, lab_filtered, on = "patient_id", how = "outer")
attrib_df = pd.merge(attrib_df, fup_filtered, on = "patient_id", how = "outer")
attrib_df = pd.merge(attrib_df, adhe_filtered, on = "patient_id", how = "outer")
cl_stats_attrib_merged = pd.merge(cl_stats, attrib_df, on = "patient_id", how = "left")

print("Exporting node attribute data to " + str(nd_attrib_file))
cl_stats_attrib_merged.to_csv(nd_attrib_file, index = False)

print("Merging attribute tables at baseline")
attrib_bl_df = pd.merge(pat_tbl, lab_filtered_bl, on = "patient_id", how = "outer")
attrib_bl_df = pd.merge(attrib_bl_df, fup_filtered_bl, on = "patient_id", how = "outer")
attrib_bl_df = pd.merge(attrib_bl_df, adhe_filtered_bl, on = "patient_id", how = "outer")
cl_stats_bl_attrib_merged = pd.merge(cl_stats_bl, attrib_bl_df, on = "patient_id", how = "left")

print("Exporting node attribute data at baseline to " + str(nd_attrib_bl_file))
cl_stats_bl_attrib_merged.to_csv(nd_attrib_bl_file, index = False)

if args.classify:
    print("Classifying clusters and exporting to " + cl_class_file)
    classify_df = classify_clusters(args.classify)
    classify_df.to_csv(cl_class_file, index = False)

#}}}

