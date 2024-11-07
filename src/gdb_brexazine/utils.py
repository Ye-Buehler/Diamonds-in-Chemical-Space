import pandas as pd
import ast
import re
import os

FILE_PATH_READ = "file_path_read"
FILE_PATH_SAVE = "file_path_save"

class Utils:

    # TODO: append all the dataframes from a folder
    def append_dfs_in_folder(self, FOLDER_PATH, COLUMN_NAME_INPUT) -> pd.DataFrame:

        iteration = 0
        row_count = 0

        for filename in os.listdir(FOLDER_PATH):
            
            iteration += 1
            if iteration == 1:
                df = pd.read_csv(FOLDER_PATH + filename, sep='\t', names=COLUMN_NAME_INPUT)
                #print("basis lenght = " + str(len(df)))
                row_count += len(df)

            if iteration != 1:
                df2 = pd.read_csv(FOLDER_PATH + filename, sep='\t',names=COLUMN_NAME_INPUT)
                #print("added lenght = " + str(len(df2)))
                df = pd.concat([df, df2], axis=0, ignore_index=True)
                row_count += len(df2)

        print("Total enteries that have been merged: \t" + str(row_count))
        return df