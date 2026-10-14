# PKIS2 Data Inspection Notes

## S4

- File: pone.0181585.s004.xlsx
- Sheet: Table 4 - PKIS2 %Inh
- Initial pandas shape: 646 rows × 413 columns
- Measurement scale appears potentially related to % inhibition based on
  the sheet name, but this still needs verification from the spreadsheet
  contents and paper Methods.

- Initial S4 shape read by pandas: 646 × 413.
- Rows 0–644 appear to contain compound records (645 rows).
- Row 645 lacks compound metadata and SMILES.
- The only non-missing value in row 645 is TEC = 13.0.
- The reason for this extra value still needs verification.