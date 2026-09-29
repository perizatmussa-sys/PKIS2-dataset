import pandas as pd


# ==========================================
# 1. FILE PATHS
# ==========================================

screening_file = "raw/pone.0181585.s004.xlsx"
kd_file = "raw/pone.0181585.s003.xlsx"


# ==========================================
# 2. FUNCTION FOR INSPECTING AN EXCEL FILE
# ==========================================

def inspect_excel(file_path):

    print("\n======================================")
    print("Inspecting file:", file_path)
    print("======================================")

    # Open Excel workbook
    excel_file = pd.ExcelFile(file_path)

    print("\nSheet names:")
    print(excel_file.sheet_names)


    # ==========================================
    # 3. INSPECT EVERY SHEET
    # ==========================================

    for sheet in excel_file.sheet_names:

        print("\n--------------------------------------")
        print("Sheet:", sheet)
        print("--------------------------------------")

        df = pd.read_excel(
            file_path,
            sheet_name=sheet
        )


        # ==========================================
        # 4. DATASET SHAPE
        # ==========================================

        print("\nShape:")
        print(df.shape)


        # ==========================================
        # 5. COLUMN NAMES
        # ==========================================

        print("\nColumns:")

        for column in df.columns:
            print(repr(column))


        # ==========================================
        # 6. FIRST 10 ROWS
        # ==========================================

        print("\nFirst 10 rows:")
        print(df.head(10).to_string())


        # ==========================================
        # 7. DATA TYPES
        # ==========================================

        print("\nData types:")
        print(df.dtypes)


        # ==========================================
        # 8. MISSING VALUES
        # ==========================================

        print("\nMissing values per column:")
        print(df.isna().sum())


        # ==========================================
        # 9. NUMBER OF UNIQUE VALUES
        # ==========================================

        print("\nUnique values per column:")
        print(df.nunique(dropna=False))


        # ==========================================
        # 10. CHECK FOR LINE BREAKS
        # ==========================================

        print("\nCells containing line breaks:")

        for column in df.columns:

            count = (
                df[column]
                .astype(str)
                .str.contains("\n", regex=False)
                .sum()
            )

            print(column, ":", count)


        # ==========================================
        # 11. DUPLICATE ROWS
        # ==========================================

        print("\nDuplicate rows:")
        print(df.duplicated().sum())


        # ==========================================
        # 12. BASIC INFORMATION
        # ==========================================

        print("\nDataFrame information:")
        df.info()


# ==========================================
# 13. INSPECT BOTH RAW FILES
# ==========================================

inspect_excel(screening_file)
inspect_excel(kd_file)