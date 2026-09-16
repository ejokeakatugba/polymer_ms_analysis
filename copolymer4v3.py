import os
import numpy as np
import pandas as pd
from tkinter import Tk, filedialog


# ============================================================
# PEG-PE COPOLYMER ANALYSIS
# VERSION 3
# ============================================================
#
# WHAT THIS VERSION DOES:
#
# 1. Reads raw m/z + intensity Excel data
#
# 2. Detects PEG homologous series:
#       +44.026215 Da
#
# 3. Detects PE homologous series:
#       +28.031300 Da
#
# 4. Keeps the DYNAMIC 5 ppm mass tolerance
#
# 5. Calculates:
#       total intensity of EACH PEG series
#       total intensity of EACH PE series
#       summed intensity of ALL PEG series
#       summed intensity of ALL PE series
#
# 6. Compares PEG series pairwise:
#
#       PEG_A       PEG_B
#         ↓ PEG       ↓ PEG
#       peak         peak
#
#       horizontal separation =
#       n × PE repeat
#
#       At least 3 corresponding members required.
#
# 7. Does reverse analysis:
#
#       PE_A        PE_B
#         ↓ PE        ↓ PE
#
#       horizontal separation =
#       n × PEG repeat
#
# 8. Saves a CLEAN Excel workbook:
#
#       Summary
#       PEG Series
#       PE Series
#       PEG Pair Candidates
#       PE Pair Candidates
#
# Candidate sheets are VISUAL rather than giant tables.
#
# ============================================================


# ============================================================
# SETTINGS
# ============================================================

PEG_REPEAT = 44.026215       # C2H4O
PE_REPEAT = 28.031300        # C2H4

PPM_TOLERANCE = 5.0

# Minimum number of peaks needed for an individual
# homologous series
MIN_SERIES_LENGTH = 4

# Minimum number of corresponding members needed before
# two series are considered a candidate parallel pair
MIN_PARALLEL_MATCHES = 3

# Search 1×, 2×, 3× ... opposite repeat
MAX_PE_MULTIPLE =8
MAX_PEG_MULTIPLE =8

# Intensity filter
BACKGROUND_INTENSITY = 500


# ============================================================
# SELECT RAW EXCEL FILE
# ============================================================

root = Tk()
root.withdraw()

file_path = filedialog.askopenfilename(
    title="Select PEG-PE raw Excel spectrum",
    filetypes=[
        ("Excel files", "*.xlsx *.xls")
    ]
)

if not file_path:
    raise SystemExit("No file selected.")


# ============================================================
# READ DATA
# ============================================================

df = pd.read_excel(file_path)

# First two columns assumed to be m/z and intensity
mz_col = df.columns[0]
intensity_col = df.columns[1]

data = df[[mz_col, intensity_col]].copy()

data.columns = [
    "m/z",
    "Intensity"
]

data["m/z"] = pd.to_numeric(
    data["m/z"],
    errors="coerce"
)

data["Intensity"] = pd.to_numeric(
    data["Intensity"],
    errors="coerce"
)

data = data.dropna()

data = data[
    data["Intensity"] > BACKGROUND_INTENSITY
]

data = data.sort_values(
    "m/z"
).reset_index(drop=True)

mz_array = data["m/z"].to_numpy(
    dtype=float
)

intensity_lookup = dict(
    zip(
        data["m/z"].round(8),
        data["Intensity"]
    )
)


# ============================================================
# DYNAMIC 5 PPM TOLERANCE
# ============================================================

def dynamic_tolerance(mz):

    return (
        mz
        * PPM_TOLERANCE
        / 1_000_000
    )


# ============================================================
# FIND HOMOLOGOUS SERIES
# ============================================================

def check_mass_dynamic(
        mz_values,
        repeat_mass):

    """
    Finds continuous homologous series.

    A second peak is accepted if:

    | observed difference - theoretical repeat |

    is within the dynamic 5 ppm tolerance
    calculated from the second peak m/z.

    This preserves the Version 2 dynamic
    mass-tolerance approach.
    """

    A = [
        float(x)
        for x in mz_values
    ]

    series_map = {}

    key_map = {}

    for i in range(len(A)):

        key = A[i]

        if key in key_map:
            key = key_map[key]

        if key not in series_map:
            series_map[key] = set()

        series_map[key].add(
            A[i]
        )

        for j in range(
                i + 1,
                len(A)):

            mz1 = A[i]
            mz2 = A[j]

            observed_difference = (
                mz2 - mz1
            )

            error = (
                observed_difference
                - repeat_mass
            )

            tolerance = (
                dynamic_tolerance(
                    mz2
                )
            )

            if abs(error) <= tolerance:

                series_map[key].add(
                    mz2
                )

                key_map[mz2] = key

    for key in series_map:

        series_map[key] = sorted(
            series_map[key]
        )

    return series_map


# ============================================================
# EXTRACT VALID SERIES
# ============================================================

def extract_series(
        series_map,
        prefix):

    series_list = []

    seen = set()

    count = 1

    for key, members in series_map.items():

        members = sorted(
            set(members)
        )

        if (
            len(members)
            < MIN_SERIES_LENGTH
        ):
            continue

        signature = tuple(
            round(x, 6)
            for x in members
        )

        if signature in seen:
            continue

        seen.add(signature)

        total_intensity = sum(

            intensity_lookup.get(
                round(mz, 8),
                0
            )

            for mz in members
        )

        series_list.append({

            "Series_ID":
                f"{prefix}_{count}",

            "Members":
                members,

            "Number_of_Peaks":
                len(members),

            "Total_Intensity":
                total_intensity
        })

        count += 1

    return series_list


# ============================================================
# FIND PEG AND PE SERIES
# ============================================================

peg_map = check_mass_dynamic(
    mz_array,
    PEG_REPEAT
)

pe_map = check_mass_dynamic(
    mz_array,
    PE_REPEAT
)

peg_series = extract_series(
    peg_map,
    "PEG"
)

pe_series = extract_series(
    pe_map,
    "PE"
)


# ============================================================
# SIMPLE SERIES TABLE
# ============================================================

def make_simple_series_table(
        series_list):

    rows = []

    for s in series_list:

        rows.append({

            "Series":
                s["Series_ID"],

            "Number of Peaks":
                s["Number_of_Peaks"],

            "Total Summed Intensity":
                s["Total_Intensity"],

            "Peak Series":
                "  →  ".join(
                    f"{mz:.4f}"
                    for mz in s["Members"]
                )
        })

    return pd.DataFrame(rows)


peg_series_df = (
    make_simple_series_table(
        peg_series
    )
)

pe_series_df = (
    make_simple_series_table(
        pe_series
    )
)


# ============================================================
# MATCH TWO PARALLEL SERIES
# ============================================================

def compare_two_series(
        series1,
        series2,
        opposite_repeat,
        max_multiple):

    """
    Example:

    PEG series vs PEG series

    Search whether corresponding members
    are separated horizontally by:

        n × PE_REPEAT

    Reverse analysis uses PE series and
    n × PEG_REPEAT.

    At least MIN_PARALLEL_MATCHES must match.
    """

    # Lower-mass series goes on LEFT
    if (
        np.mean(
            series1["Members"]
        )
        <=
        np.mean(
            series2["Members"]
        )
    ):

        left_series = series1
        right_series = series2

    else:

        left_series = series2
        right_series = series1

    best_result = None

    # --------------------------------------------------------
    # Try 1×, 2×, 3× ... opposite repeat
    # --------------------------------------------------------

    for n in range(
            1,
            max_multiple + 1):

        expected_difference = (
            n
            * opposite_repeat
        )

        possible_matches = []

        for left_mz in (
                left_series["Members"]):

            for right_mz in (
                    right_series["Members"]):

                if (
                    right_mz
                    <=
                    left_mz
                ):
                    continue

                observed_difference = (
                    right_mz
                    - left_mz
                )

                error = (
                    observed_difference
                    - expected_difference
                )

                tolerance = (
                    dynamic_tolerance(
                        right_mz
                    )
                )

                if (
                    abs(error)
                    <=
                    tolerance
                ):

                    possible_matches.append({

                        "Left_mz":
                            left_mz,

                        "Right_mz":
                            right_mz,

                        "Observed_Difference":
                            observed_difference,

                        "Expected_Difference":
                            expected_difference,

                        "Error_Da":
                            error,

                        "Tolerance_Da":
                            tolerance
                    })

        if not possible_matches:
            continue

        # ----------------------------------------------------
        # One peak cannot be matched to several peaks.
        #
        # Sort by smallest mass error and retain one-to-one
        # matches.
        # ----------------------------------------------------

        possible_matches = sorted(

            possible_matches,

            key=lambda x:
                abs(
                    x["Error_Da"]
                )
        )

        used_left = set()
        used_right = set()

        accepted_matches = []

        for match in possible_matches:

            left_key = round(
                match["Left_mz"],
                8
            )

            right_key = round(
                match["Right_mz"],
                8
            )

            if left_key in used_left:
                continue

            if right_key in used_right:
                continue

            used_left.add(
                left_key
            )

            used_right.add(
                right_key
            )

            accepted_matches.append(
                match
            )

        accepted_matches = sorted(

            accepted_matches,

            key=lambda x:
                x["Left_mz"]
        )

        # ----------------------------------------------------
        # REQUIRE AT LEAST 3 CORRESPONDING MEMBERS
        # ----------------------------------------------------

        if (
            len(accepted_matches)
            <
            MIN_PARALLEL_MATCHES
        ):
            continue

        mean_error = np.mean([

            abs(
                x["Error_Da"]
            )

            for x
            in accepted_matches
        ])

        result = {

            "Left_Series":
                left_series,

            "Right_Series":
                right_series,

            "Multiple":
                n,

            "Expected_Difference":
                expected_difference,

            "Matches":
                accepted_matches,

            "Number_of_Matches":
                len(
                    accepted_matches
                ),

            "Mean_Absolute_Error":
                mean_error
        }

        # ----------------------------------------------------
        # If several n values work:
        #
        # 1. prefer more corresponding peaks
        # 2. then prefer smaller mass error
        # ----------------------------------------------------

        if best_result is None:

            best_result = result

        else:

            if (
                result[
                    "Number_of_Matches"
                ]
                >
                best_result[
                    "Number_of_Matches"
                ]
            ):

                best_result = result

            elif (
                result[
                    "Number_of_Matches"
                ]
                ==
                best_result[
                    "Number_of_Matches"
                ]
            ):

                if (
                    result[
                        "Mean_Absolute_Error"
                    ]
                    <
                    best_result[
                        "Mean_Absolute_Error"
                    ]
                ):

                    best_result = result

    return best_result


# ============================================================
# COMPARE EVERY POSSIBLE SERIES PAIR
# ============================================================

def find_parallel_candidates(
        series_list,
        opposite_repeat,
        max_multiple):

    candidates = []

    for i in range(
            len(series_list)):

        for j in range(
                i + 1,
                len(series_list)):

            result = (
                compare_two_series(

                    series_list[i],

                    series_list[j],

                    opposite_repeat,

                    max_multiple
                )
            )

            if result is not None:

                candidates.append(
                    result
                )

    # --------------------------------------------------------
    # Sort:
    # candidates with most matching members first
    # --------------------------------------------------------

    candidates = sorted(

        candidates,

        key=lambda x:
            (
                -x[
                    "Number_of_Matches"
                ],

                x[
                    "Mean_Absolute_Error"
                ]
            )
    )

    return candidates


# ============================================================
# PEG PAIR CANDIDATES
#
# PEG series internally differ vertically by PEG.
#
# Two PEG series must be horizontally separated by PE.
# ============================================================

peg_pair_candidates = (
    find_parallel_candidates(

        peg_series,

        PE_REPEAT,

        MAX_PE_MULTIPLE
    )
)


# ============================================================
# PE PAIR CANDIDATES
#
# PE series internally differ vertically by PE.
#
# Two PE series must be horizontally separated by PEG.
# ============================================================

pe_pair_candidates = (
    find_parallel_candidates(

        pe_series,

        PEG_REPEAT,

        MAX_PEG_MULTIPLE
    )
)


# ============================================================
# TOTAL SERIES INTENSITIES
# ============================================================

total_peg_series_intensity = sum(

    s["Total_Intensity"]

    for s
    in peg_series
)

total_pe_series_intensity = sum(

    s["Total_Intensity"]

    for s
    in pe_series
)


# ============================================================
# SUMMARY
# ============================================================

summary_df = pd.DataFrame({

    "Metric": [

        "Dynamic mass tolerance",

        "PEG repeat mass (C2H4O)",

        "PE repeat mass (C2H4)",

        "Minimum peaks per series",

        "Minimum corresponding peaks for pair candidate",

        "PEG series detected",

        "PE series detected",

        "Total summed intensity of PEG series",

        "Total summed intensity of PE series",

        "PEG pair candidates",

        "PE pair candidates"
    ],

    "Value": [

        f"{PPM_TOLERANCE} ppm",

        PEG_REPEAT,

        PE_REPEAT,

        MIN_SERIES_LENGTH,

        MIN_PARALLEL_MATCHES,

        len(peg_series),

        len(pe_series),

        total_peg_series_intensity,

        total_pe_series_intensity,

        len(
            peg_pair_candidates
        ),

        len(
            pe_pair_candidates
        )
    ]
})


# ============================================================
# OUTPUT FILE
# ============================================================

folder = os.path.dirname(
    file_path
)

base_name = os.path.splitext(
    os.path.basename(
        file_path
    )
)[0]

output_path = os.path.join(

    folder,

    base_name
    + "_VERSION3_PARALLEL_COPOLYMER.xlsx"
)


# ============================================================
# SAVE WORKBOOK
# ============================================================

with pd.ExcelWriter(
        output_path,
        engine="xlsxwriter"
) as writer:

    workbook = writer.book


    # ========================================================
    # FORMATS
    # ========================================================

    title_format = (
        workbook.add_format({

            "bold": True,
            "font_size": 16,
            "align": "center",
            "valign": "vcenter",
            "bg_color": "#D9EAF7",
            "border": 1
        })
    )

    candidate_format = (
        workbook.add_format({

            "bold": True,
            "font_size": 14,
            "font_color": "#C00000"
        })
    )

    series_name_format = (
        workbook.add_format({

            "bold": True,
            "font_size": 12,
            "align": "center"
        })
    )

    intensity_format = (
        workbook.add_format({

            "italic": True,
            "align": "center"
        })
    )

    peak_format = (
        workbook.add_format({

            "font_size": 12,
            "align": "center",
            "valign": "vcenter",
            "border": 1
        })
    )

    arrow_format = (
        workbook.add_format({

            "bold": True,
            "font_size": 11,
            "align": "center",
            "font_color": "#ED7D31"
        })
    )

    vertical_format = (
        workbook.add_format({

            "bold": True,
            "align": "center",
            "font_color": "#1F4E78"
        })
    )

    info_format = (
        workbook.add_format({

            "bold": True,
            "font_size": 11
        })
    )

    header_format = (
        workbook.add_format({

            "bold": True,
            "bg_color": "#D9EAD3",
            "border": 1,
            "align": "center"
        })
    )


    # ========================================================
    # SUMMARY
    # ========================================================

    summary_df.to_excel(

        writer,

        sheet_name="Summary",

        index=False
    )

    summary_ws = (
        writer.sheets[
            "Summary"
        ]
    )

    summary_ws.set_column(
        "A:A",
        48
    )

    summary_ws.set_column(
        "B:B",
        24
    )

    summary_ws.set_row(
        0,
        24
    )


    # ========================================================
    # PEG SERIES
    # ========================================================

    peg_series_df.to_excel(

        writer,

        sheet_name="PEG Series",

        index=False
    )

    peg_ws = (
        writer.sheets[
            "PEG Series"
        ]
    )

    peg_ws.set_column(
        "A:A",
        12
    )

    peg_ws.set_column(
        "B:B",
        18
    )

    peg_ws.set_column(
        "C:C",
        24
    )

    peg_ws.set_column(
        "D:D",
        80
    )


    # ========================================================
    # PE SERIES
    # ========================================================

    pe_series_df.to_excel(

        writer,

        sheet_name="PE Series",

        index=False
    )

    pe_ws = (
        writer.sheets[
            "PE Series"
        ]
    )

    pe_ws.set_column(
        "A:A",
        12
    )

    pe_ws.set_column(
        "B:B",
        18
    )

    pe_ws.set_column(
        "C:C",
        24
    )

    pe_ws.set_column(
        "D:D",
        80
    )


    # ========================================================
    # FUNCTION TO DRAW CANDIDATES
    # ========================================================

    def write_candidate_sheet(
            worksheet,
            candidates,
            candidate_type):

        worksheet.set_column(
            "A:A",
            23
        )

        worksheet.set_column(
            "B:B",
            30
        )

        worksheet.set_column(
            "C:C",
            23
        )

        worksheet.set_column(
            "D:D",
            3
        )

        worksheet.set_column(
            "E:E",
            30
        )

        if candidate_type == "PEG":

            vertical_name = (
                "C₂H₄O"
            )

            vertical_mass = (
                PEG_REPEAT
            )

            horizontal_name = (
                "C₂H₄"
            )

        else:

            vertical_name = (
                "C₂H₄"
            )

            vertical_mass = (
                PE_REPEAT
            )

            horizontal_name = (
                "C₂H₄O"
            )


        # ----------------------------------------------------
        # NO CANDIDATES
        # ----------------------------------------------------

        if len(candidates) == 0:

            worksheet.write(
                0,
                0,
                "No qualifying parallel-series candidates found."
            )

            return


        row = 0


        # ----------------------------------------------------
        # EACH CANDIDATE GETS ITS OWN VISUAL BLOCK
        # ----------------------------------------------------

        for candidate_number, candidate in enumerate(
                candidates,
                start=1):

            left = (
                candidate[
                    "Left_Series"
                ]
            )

            right = (
                candidate[
                    "Right_Series"
                ]
            )

            n = (
                candidate[
                    "Multiple"
                ]
            )

            expected_difference = (
                candidate[
                    "Expected_Difference"
                ]
            )

            matches = (
                candidate[
                    "Matches"
                ]
            )


            # =================================================
            # CANDIDATE TITLE
            # =================================================

            worksheet.merge_range(

                row,
                0,
                row,
                2,

                (
                    f"{candidate_type} PAIR "
                    f"CANDIDATE {candidate_number}"
                ),

                candidate_format
            )

            row += 1


            # =================================================
            # RELATIONSHIP INFORMATION
            # =================================================

            worksheet.merge_range(

                row,
                0,
                row,
                2,

                (
                    f"{left['Series_ID']}  ↔  "
                    f"{right['Series_ID']}"
                ),

                title_format
            )

            row += 1


            worksheet.merge_range(

                row,
                0,
                row,
                2,

                (
                    f"Horizontal relationship: "
                    f"{n} × {horizontal_name} "
                    f"= {expected_difference:.6f} Da"
                ),

                info_format
            )

            row += 1


            worksheet.merge_range(

                row,
                0,
                row,
                2,

                (
                    f"Matched corresponding members: "
                    f"{len(matches)}"
                ),

                info_format
            )

            row += 2


            # =================================================
            # SERIES HEADERS
            # =================================================

            worksheet.write(

                row,
                0,
                left[
                    "Series_ID"
                ],

                series_name_format
            )

            worksheet.write(

                row,
                2,
                right[
                    "Series_ID"
                ],

                series_name_format
            )

            row += 1


            # =================================================
            # SERIES TOTAL INTENSITIES
            # =================================================

            worksheet.write(

                row,
                0,

                (
                    "Total intensity = "
                    f"{left['Total_Intensity']:.4f}"
                ),

                intensity_format
            )

            worksheet.write(

                row,
                2,

                (
                    "Total intensity = "
                    f"{right['Total_Intensity']:.4f}"
                ),

                intensity_format
            )

            row += 2


            # =================================================
            # PEAKS
            # =================================================

            for match_number, match in enumerate(
                    matches):

                left_mz = (
                    match[
                        "Left_mz"
                    ]
                )

                right_mz = (
                    match[
                        "Right_mz"
                    ]
                )

                observed = (
                    match[
                        "Observed_Difference"
                    ]
                )


                # ---------------------------------------------
                # LEFT PEAK
                # ---------------------------------------------

                worksheet.write(

                    row,
                    0,

                    f"{left_mz:.4f}",

                    peak_format
                )


                # ---------------------------------------------
                # HORIZONTAL ARROW
                # ---------------------------------------------

                worksheet.write(

                    row,
                    1,

                    (
                        f"── {horizontal_name} →\n"
                        f"Δ = {observed:.4f}"
                    ),

                    arrow_format
                )


                # ---------------------------------------------
                # RIGHT PEAK
                # ---------------------------------------------

                worksheet.write(

                    row,
                    2,

                    f"{right_mz:.4f}",

                    peak_format
                )


                worksheet.set_row(
                    row,
                    32
                )


                # ---------------------------------------------
                # VERTICAL ARROWS BETWEEN MEMBERS
                # ---------------------------------------------

                if (
                    match_number
                    <
                    len(matches) - 1
                ):

                    row += 1

                    worksheet.write(

                        row,
                        0,

                        (
                            f"↓  +{vertical_name}\n"
                            f"{vertical_mass:.4f} Da"
                        ),

                        vertical_format
                    )

                    worksheet.write(

                        row,
                        2,

                        (
                            f"↓  +{vertical_name}\n"
                            f"{vertical_mass:.4f} Da"
                        ),

                        vertical_format
                    )

                    worksheet.set_row(
                        row,
                        28
                    )

                row += 1


            # =================================================
            # MASS ERROR INFO
            #
            # Kept simple. Dynamic tolerance still used in
            # background to ACCEPT/REJECT relationships.
            # =================================================

            row += 1

            worksheet.merge_range(

                row,
                0,
                row,
                2,

                (
                    f"Dynamic matching criterion: "
                    f"{PPM_TOLERANCE} ppm"
                ),

                info_format
            )

            row += 3


    # ========================================================
    # PEG PAIR CANDIDATES
    # ========================================================

    peg_pair_ws = (
        workbook.add_worksheet(
            "PEG Pair Candidates"
        )
    )

    writer.sheets[
        "PEG Pair Candidates"
    ] = peg_pair_ws

    write_candidate_sheet(

        peg_pair_ws,

        peg_pair_candidates,

        "PEG"
    )


    # ========================================================
    # PE PAIR CANDIDATES
    # ========================================================

    pe_pair_ws = (
        workbook.add_worksheet(
            "PE Pair Candidates"
        )
    )

    writer.sheets[
        "PE Pair Candidates"
    ] = pe_pair_ws

    write_candidate_sheet(

        pe_pair_ws,

        pe_pair_candidates,

        "PE"
    )


# ============================================================
# CONSOLE SUMMARY
# ============================================================

print("\n")
print("=" * 60)
print("VERSION 3 ANALYSIS COMPLETE")
print("=" * 60)

print(
    f"\nPEG series detected: "
    f"{len(peg_series)}"
)

print(
    f"PE series detected: "
    f"{len(pe_series)}"
)

print(
    "\nTotal PEG-series intensity: "
    f"{total_peg_series_intensity:.4f}"
)

print(
    "Total PE-series intensity: "
    f"{total_pe_series_intensity:.4f}"
)

print(
    "\nPEG pair candidates: "
    f"{len(peg_pair_candidates)}"
)

print(
    "PE pair candidates: "
    f"{len(pe_pair_candidates)}"
)

print(
    f"\nSaved workbook:\n"
    f"{output_path}"
)

print("\n")