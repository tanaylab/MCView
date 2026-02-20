#!/usr/bin/env python3
"""
Metacell Import Script

This script converts data from the metacell1 format to the newer metacells format.
It processes single-cell data to generate properly formatted metacell data for further analysis.
Optionally, it can:
- Import custom x,y coordinates for the metacells from a CSV file
- Import a graph structure from a CSV file to set as obs_outgoing_weights

Usage:
    python mc1_to_mc2.py --input cells.h5ad --output metacells.h5ad [--coords coords.csv] [--graph graph.csv] [--seed 60427]
"""

import sys
sys.path = ['/home/ofirr/miniforge3/envs/amosbase8/bin',
     '/home/ofirr/miniforge3/envs/amosbase8/lib/python312.zip',
    '/home/ofirr/miniforge3/envs/amosbase8/lib/python3.12',
    '/home/ofirr/miniforge3/envs/amosbase8/lib/python3.12/lib-dynload',
    '/home/ofirr/miniforge3/envs/amosbase8/lib/python3.12/site-packages']

import argparse
import sys
import logging
import metacells as mc
import anndata as ad
import pandas as pd
import os
import numpy as np
from scipy import sparse


def setup_logging():
    """Configure logging for the script."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Import and convert metacell1 data to the newer metacells format')
    
    parser.add_argument('--input', '-i', required=True, 
                        help='Input h5ad file containing metacell1 format data')
    
    parser.add_argument('--output', '-o', required=True,
                        help='Output h5ad file for storing converted metacells data')
    
    parser.add_argument('--coords', '-c', 
                        help='Optional CSV file with metacell,x,y coordinates to override the default layout')
    
    parser.add_argument('--graph', '-g',
                        help='Optional CSV file with from,to,weight defining the graph structure')
    
    parser.add_argument('--seed', '-s', type=int, default=60427,
                        help='Random seed for reproducibility (default: 60427)')
    
    return parser.parse_args()


def import_coordinates(coords_file):
    """Import metacell coordinates from a CSV file.
    
    Args:
        coords_file: Path to a CSV file with columns: metacell, x, y
        
    Returns:
        DataFrame with metacell coordinates
    """
    if not os.path.exists(coords_file):
        raise FileNotFoundError(f"Coordinates file not found: {coords_file}")
    
    try:
        coords_df = pd.read_csv(coords_file)
        required_columns = ['metacell', 'x', 'y']
        
        for col in required_columns:
            if col not in coords_df.columns:
                raise ValueError(f"Coordinates file must contain column: {col}")
        
        # Set metacell column as index
        coords_df = coords_df.set_index('metacell')
        return coords_df
    
    except Exception as e:
        raise ValueError(f"Error parsing coordinates file: {e}")


def import_graph(graph_file, metacell_ids):
    """Import graph structure from a CSV file and convert to a sparse matrix.
    
    Args:
        graph_file: Path to a CSV file with columns: from, to, weight
        metacell_ids: List of metacell IDs to use as indices for the sparse matrix
        
    Returns:
        A scipy sparse matrix representing the graph
    """
    if not os.path.exists(graph_file):
        raise FileNotFoundError(f"Graph file not found: {graph_file}")
    
    try:
        graph_df = pd.read_csv(graph_file)
        required_columns = ['from', 'to', 'weight']
        
        for col in required_columns:
            if col not in graph_df.columns:
                raise ValueError(f"Graph file must contain column: {col}")
        
        # Create a mapping of metacell IDs to indices
        metacell_to_idx = {mc_id: idx for idx, mc_id in enumerate(metacell_ids)}
        
        # Extract row, col, and data for sparse matrix
        rows = []
        cols = []
        data = []
        
        for _, row in graph_df.iterrows():
            from_metacell = row['from']
            to_metacell = row['to']
            weight = row['weight']
            
            if from_metacell in metacell_to_idx and to_metacell in metacell_to_idx:
                rows.append(metacell_to_idx[from_metacell])
                cols.append(metacell_to_idx[to_metacell])
                data.append(weight)
            
        # Create the sparse matrix
        n = len(metacell_ids)
        sparse_matrix = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))
        
        return sparse_matrix
    
    except Exception as e:
        raise ValueError(f"Error parsing graph file: {e}")


def convert_metacells(input_file, output_file, coords_file=None, graph_file=None, random_seed=60427):
    """Convert data from metacell1 format to metacells format.
    
    Args:
        input_file: Path to the input h5ad file in metacell1 format
        output_file: Path to the output h5ad file in metacells format
        coords_file: Optional path to a CSV file with metacell coordinates
        graph_file: Optional path to a CSV file with graph structure
        random_seed: Integer seed for reproducibility
    """
    logging.info(f"Reading metacell1 data from {input_file}")
    cdata = ad.read_h5ad(input_file)
    
    logging.info("Collecting metacells in new format")
    mdata = mc.pl.collect_metacells(cdata, name="metacells", random_seed=random_seed)
    
    logging.info("Computing data for mcview visualization")
    mc.pl.compute_for_mcview(adata=cdata, gdata=mdata, random_seed=random_seed)
    
    logging.info("Converting metacell naming conventions")
    metacell_mapping = cdata.obs[['metacell_name', 'metacell_orig']].drop_duplicates().set_index('metacell_name')['metacell_orig'].to_dict()
    
    cdata.obs.drop('metacell_name', axis=1, inplace=True)
    cdata.obs.rename(columns={'metacell_orig': 'metacell_name'}, inplace=True)
    
    mdata.obs['metacell_name'] = [metacell_mapping[i] for i in mdata.obs.index]
    mdata.obs.index = mdata.obs.metacell_name
    
    # Import custom coordinates if provided
    if coords_file:
        logging.info(f"Importing custom coordinates from {coords_file}")
        try:
            coords_df = import_coordinates(coords_file)
            
            # Check if all metacells have coordinates
            missing_metacells = set(mdata.obs.index) - set(coords_df.index)
            if missing_metacells:
                logging.warning(f"Missing coordinates for {len(missing_metacells)} metacells. "
                               f"First few missing: {list(missing_metacells)[:5]}")
            
            # Update x and y coordinates in mdata.obs
            mdata.obs.loc[coords_df.index, 'x'] = coords_df['x']
            mdata.obs.loc[coords_df.index, 'y'] = coords_df['y']
            
            logging.info(f"Updated coordinates for {len(coords_df)} metacells")
        except Exception as e:
            logging.error(f"Failed to import coordinates: {e}")
            logging.info("Continuing with default coordinates")
    
    # Import graph structure if provided
    if graph_file:
        logging.info(f"Importing graph structure from {graph_file}")
        try:
            # Get ordered list of metacell IDs to use as indices for the sparse matrix
            metacell_ids = list(mdata.obs.index)
            
            # Import graph as sparse matrix
            graph_matrix = import_graph(graph_file, metacell_ids)
            
            # Set the graph as obs_outgoing_weights
            if 'obsp' not in mdata:
                mdata.obsp = {}
                
            mdata.obsp['obs_outgoing_weights'] = graph_matrix
            
            logging.info(f"Set graph structure with {graph_matrix.count_nonzero()} edges")
        except Exception as e:
            logging.error(f"Failed to import graph: {e}")
            logging.info("Continuing without custom graph structure")
    
    logging.info(f"Writing converted metacells data to {output_file}")
    mdata.write_h5ad(output_file)
    
    logging.info("Conversion complete")


def main():
    """Main entry point of the script."""
    setup_logging()
    args = parse_arguments()
    
    try:
        convert_metacells(args.input, args.output, args.coords, args.graph, args.seed)
    except Exception as e:
        logging.error(f"Error converting metacells data: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())