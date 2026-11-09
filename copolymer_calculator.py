# ============================================================
# GENERAL COPOLYMER EXACT-MASS ASSIGNMENT CALCULATOR
# Dash application
# ============================================================

from dash import Dash, html, dcc, Input, Output, State, dash_table
import pandas as pd


# ============================================================
# CREATE APP
# ============================================================

app = Dash(__name__)
app.title = "Copolymer Mass Calculator"


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def ppm_error(observed_mz, calculated_mz):
    """Calculate signed mass error in ppm."""
    return ((observed_mz - calculated_mz) / calculated_mz) * 1_000_000


def parse_mass_list(text):
    """
    Convert text such as:

    C2H6O,46.041865
    OH,17.002740

    into:
    [
        ("C2H6O", 46.041865),
        ("OH", 17.002740)
    ]
    """

    results = []

    if not text:
        return results

    for line in text.splitlines():

        line = line.strip()

        if not line:
            continue

        try:
            name, mass = line.rsplit(",", 1)

            name = name.strip()
            mass = float(mass.strip())

            results.append((name, mass))

        except ValueError:
            continue

    return results


def search_assignments(
    observed_mz,
    repeat_a_name,
    repeat_a_mass,
    repeat_b_name,
    repeat_b_mass,
    end_groups,
    adducts,
    ppm_tolerance,
    max_a,
    max_b
):

    results = []

    # Is this a copolymer?
    use_repeat_b = (
        repeat_b_name is not None
        and str(repeat_b_name).strip() != ""
        and repeat_b_mass is not None
    )

    if not use_repeat_b:
        repeat_b_name = ""
        repeat_b_mass = 0
        max_b = 0

    # --------------------------------------------------------
    # Search every possible repeat combination
    # --------------------------------------------------------

    for a_count in range(max_a + 1):

        for b_count in range(max_b + 1):

            # Must contain at least one polymer repeat
            if a_count == 0 and b_count == 0:
                continue

            repeat_mass = (
                a_count * repeat_a_mass
                + b_count * repeat_b_mass
            )

            for end_name, end_mass in end_groups:

                for adduct_name, adduct_mass in adducts:

                    calculated_mz = (
                        repeat_mass
                        + end_mass
                        + adduct_mass
                    )

                    if calculated_mz <= 0:
                        continue

                    error = ppm_error(
                        observed_mz,
                        calculated_mz
                    )

                    if abs(error) <= ppm_tolerance:

                        # ------------------------------------
                        # Create assignment label
                        # ------------------------------------

                        assignment = (
                            f"{repeat_a_name}{a_count}"
                        )

                        if use_repeat_b:
                            assignment += (
                                f"{repeat_b_name}{b_count}"
                            )

                        # Classify fragment
                        if use_repeat_b:

                            if a_count > 0 and b_count > 0:
                                fragment_type = "Mixed copolymer"

                            elif a_count > 0:
                                fragment_type = (
                                    f"{repeat_a_name}-only"
                                )

                            else:
                                fragment_type = (
                                    f"{repeat_b_name}-only"
                                )

                        else:
                            fragment_type = (
                                f"{repeat_a_name} fragment"
                            )

                        results.append({

                            "Assignment":
                                assignment,

                            "Fragment type":
                                fragment_type,

                            f"{repeat_a_name} count":
                                a_count,

                            (
                                f"{repeat_b_name} count"
                                if use_repeat_b
                                else "Repeat B count"
                            ):
                                b_count,

                            "End group":
                                end_name,

                            "Adduct / transformation":
                                adduct_name,

                            "Calculated m/z":
                                round(calculated_mz, 6),

                            "Observed m/z":
                                round(observed_mz, 6),

                            "ppm error":
                                round(error, 2),

                            "Absolute ppm error":
                                abs(error)
                        })

    # --------------------------------------------------------
    # Best mass match first
    # --------------------------------------------------------

    results.sort(
        key=lambda x: x["Absolute ppm error"]
    )

    # Remove internal sorting value
    for result in results:
        result.pop("Absolute ppm error", None)

    return results


# ============================================================
# DEFAULT EXAMPLE VALUES
#
# These are ONLY examples.
# User can replace them with any copolymer.
# ============================================================

DEFAULT_END_GROUPS = """C6H8O2,112.052430
C2H4O2,60.021130
C8H12,108.093900
C2H6O (combined termini),46.041865
H2O / OH+H,18.010565
OH,17.002740"""


DEFAULT_ADDUCTS = """[M+H]+,1.007276
M radical cation,-0.000549
[M+NH4]+,18.033823
[M+H+O]+,17.002191
[M+O2]+,31.989281
[M+NO]+,29.997440
[M+NO2]+,45.992354
[M+O+H-H2O]+,-1.008374
[M-H2O],-18.010565
[M-CO],-27.994915
[M-CO2],-43.989829
[M+H+H2O]+,19.017841
[M+H+2H2O]+,37.028406
[M+Na]+,22.989218
[M+K]+,38.963158"""


# ============================================================
# STYLE
# ============================================================

input_style = {
    "width": "100%",
    "padding": "10px",
    "marginBottom": "12px",
    "border": "1px solid #bbb",
    "borderRadius": "6px",
    "boxSizing": "border-box"
}


section_style = {
    "backgroundColor": "#f7f7f7",
    "padding": "20px",
    "borderRadius": "10px",
    "marginBottom": "20px"
}


label_style = {
    "fontWeight": "bold",
    "marginBottom": "5px",
    "display": "block"
}


# ============================================================
# APP LAYOUT
# ============================================================

app.layout = html.Div(

    style={
        "maxWidth": "1200px",
        "margin": "30px auto",
        "fontFamily": "Arial, sans-serif",
        "padding": "20px"
    },

    children=[

        html.H1(
            "Copolymer Exact-Mass Assignment Calculator"
        ),

        html.P(
            """
            Search possible polymer/copolymer compositions from
            an experimental m/z using user-defined repeating
            units, end groups, adducts and source transformations.
            """
        ),

        # ====================================================
        # OBSERVED MASS
        # ====================================================

        html.Div(
            style=section_style,
            children=[

                html.H3("Experimental Peak"),

                html.Label(
                    "Observed m/z",
                    style=label_style
                ),

                dcc.Input(
                    id="observed-mz",
                    type="number",
                    value=815.7691,
                    step=0.0001,
                    style=input_style
                ),

                html.Label(
                    "Mass tolerance (ppm)",
                    style=label_style
                ),

                dcc.Input(
                    id="ppm-tolerance",
                    type="number",
                    value=5,
                    min=0.1,
                    step=0.1,
                    style=input_style
                )
            ]
        ),

        # ====================================================
        # REPEATING UNITS
        # ====================================================

        html.Div(
            style=section_style,
            children=[

                html.H3("Repeating Units"),

                html.P(
                    """
                    Repeat Unit B is optional. Leave it blank
                    when analyzing a homopolymer.
                    """
                ),

                html.Div(
                    style={
                        "display": "grid",
                        "gridTemplateColumns": "1fr 1fr",
                        "gap": "20px"
                    },

                    children=[

                        html.Div([

                            html.H4("Repeat Unit A"),

                            html.Label(
                                "Name",
                                style=label_style
                            ),

                            dcc.Input(
                                id="repeat-a-name",
                                value="PE",
                                type="text",
                                placeholder="e.g. PE",
                                style=input_style
                            ),

                            html.Label(
                                "Exact repeat mass (Da)",
                                style=label_style
                            ),

                            dcc.Input(
                                id="repeat-a-mass",
                                value=28.031300,
                                type="number",
                                step=0.000001,
                                style=input_style
                            ),

                            html.Label(
                                "Maximum repeat count",
                                style=label_style
                            ),

                            dcc.Input(
                                id="max-a",
                                value=50,
                                type="number",
                                min=1,
                                step=1,
                                style=input_style
                            )

                        ]),

                        html.Div([

                            html.H4(
                                "Repeat Unit B (optional)"
                            ),

                            html.Label(
                                "Name",
                                style=label_style
                            ),

                            dcc.Input(
                                id="repeat-b-name",
                                value="EO",
                                type="text",
                                placeholder="e.g. EO",
                                style=input_style
                            ),

                            html.Label(
                                "Exact repeat mass (Da)",
                                style=label_style
                            ),

                            dcc.Input(
                                id="repeat-b-mass",
                                value=44.026215,
                                type="number",
                                step=0.000001,
                                style=input_style
                            ),

                            html.Label(
                                "Maximum repeat count",
                                style=label_style
                            ),

                            dcc.Input(
                                id="max-b",
                                value=50,
                                type="number",
                                min=0,
                                step=1,
                                style=input_style
                            )

                        ])
                    ]
                )
            ]
        ),

        # ====================================================
        # END GROUPS
        # ====================================================

        html.Div(
            style=section_style,
            children=[

                html.H3("Possible End Groups"),

                html.P(
                    "Enter one end group per line as:"
                ),

                html.Code(
                    "Name,ExactMass"
                ),

                html.Br(),
                html.Br(),

                dcc.Textarea(
                    id="end-groups",
                    value=DEFAULT_END_GROUPS,
                    style={
                        "width": "100%",
                        "height": "180px",
                        "fontFamily": "monospace",
                        "padding": "10px",
                        "boxSizing": "border-box"
                    }
                )
            ]
        ),

        # ====================================================
        # ADDUCTS
        # ====================================================

        html.Div(
            style=section_style,
            children=[

                html.H3(
                    "Adducts / Source Transformations"
                ),

                html.P([
                    """
                    Enter one possibility per line as
                    Name,MassShift. Use a negative mass for
                    a neutral loss.
                    """
                ]),

                html.Br(),

                dcc.Textarea(
                    id="adducts",
                    value=DEFAULT_ADDUCTS,
                    style={
                        "width": "100%",
                        "height": "300px",
                        "fontFamily": "monospace",
                        "padding": "10px",
                        "boxSizing": "border-box"
                    }
                )
            ]
        ),

        # ====================================================
        # CALCULATE
        # ====================================================

        html.Button(
            "Calculate Assignments",
            id="calculate-button",
            n_clicks=0,
            style={
                "padding": "14px 30px",
                "fontSize": "17px",
                "fontWeight": "bold",
                "cursor": "pointer",
                "marginBottom": "25px"
            }
        ),

        # ====================================================
        # STATUS MESSAGE
        # ====================================================

        html.Div(
            id="assignment-message",
            style={
                "fontSize": "17px",
                "fontWeight": "bold",
                "marginBottom": "20px"
            }
        ),

        # ====================================================
        # RESULTS
        # ====================================================

        html.Div(
            id="assignment-results"
        )
    ]
)


# ============================================================
# CALLBACK
# ============================================================

@app.callback(

    Output(
        "assignment-message",
        "children"
    ),

    Output(
        "assignment-results",
        "children"
    ),

    Input(
        "calculate-button",
        "n_clicks"
    ),

    State(
        "observed-mz",
        "value"
    ),

    State(
        "ppm-tolerance",
        "value"
    ),

    State(
        "repeat-a-name",
        "value"
    ),

    State(
        "repeat-a-mass",
        "value"
    ),

    State(
        "repeat-b-name",
        "value"
    ),

    State(
        "repeat-b-mass",
        "value"
    ),

    State(
        "max-a",
        "value"
    ),

    State(
        "max-b",
        "value"
    ),

    State(
        "end-groups",
        "value"
    ),

    State(
        "adducts",
        "value"
    )
)

def calculate_assignments(
    n_clicks,
    observed_mz,
    ppm_tolerance,
    repeat_a_name,
    repeat_a_mass,
    repeat_b_name,
    repeat_b_mass,
    max_a,
    max_b,
    end_group_text,
    adduct_text
):

    # --------------------------------------------------------
    # Wait until button is clicked
    # --------------------------------------------------------

    if n_clicks == 0:
        return (
            "Enter the chemistry and click Calculate Assignments.",
            ""
        )

    # --------------------------------------------------------
    # Validate basic inputs
    # --------------------------------------------------------

    if observed_mz is None:
        return (
            "Please enter an observed m/z.",
            ""
        )

    if not repeat_a_name or repeat_a_mass is None:
        return (
            "Repeat Unit A requires a name and exact mass.",
            ""
        )

    if ppm_tolerance is None:
        ppm_tolerance = 5.0

    if max_a is None:
        max_a = 50

    if max_b is None:
        max_b = 50

    # --------------------------------------------------------
    # Parse chemistry lists
    # --------------------------------------------------------

    end_groups = parse_mass_list(
        end_group_text
    )

    adducts = parse_mass_list(
        adduct_text
    )

    if not end_groups:
        return (
            "Please enter at least one valid end group.",
            ""
        )

    if not adducts:
        return (
            """
            Please enter at least one valid adduct or
            source transformation.
            """,
            ""
        )

    # --------------------------------------------------------
    # SEARCH
    # --------------------------------------------------------

    results = search_assignments(

        observed_mz=float(observed_mz),

        repeat_a_name=str(
            repeat_a_name
        ).strip(),

        repeat_a_mass=float(
            repeat_a_mass
        ),

        repeat_b_name=(
            str(repeat_b_name).strip()
            if repeat_b_name
            else None
        ),

        repeat_b_mass=(
            float(repeat_b_mass)
            if repeat_b_mass is not None
            else None
        ),

        end_groups=end_groups,

        adducts=adducts,

        ppm_tolerance=float(
            ppm_tolerance
        ),

        max_a=int(max_a),

        max_b=int(max_b)
    )

    # --------------------------------------------------------
    # NO MATCH
    # --------------------------------------------------------

    if not results:

        message = (
            f"No valid assignment for m/z "
            f"{observed_mz:.4f} within "
            f"±{ppm_tolerance} ppm. "
            f"Exact composition unresolved."
        )

        return message, ""

    # --------------------------------------------------------
    # CREATE DATAFRAME
    # --------------------------------------------------------

    df = pd.DataFrame(results)

    # Limit display to best 100 candidates
    df = df.head(100)

    best = df.iloc[0]

    message = html.Div([

        html.Div(
            f"{len(results)} candidate assignment(s) "
            f"found within ±{ppm_tolerance} ppm."
        ),

        html.Div(
            [
                "Best mass match: ",
                html.B(
                    best["Assignment"]
                ),
                " | ",
                best["End group"],
                " | ",
                best[
                    "Adduct / transformation"
                ],
                " | ",
                f"{best['ppm error']:+.2f} ppm"
            ],
            style={
                "marginTop": "8px"
            }
        ),

        html.Div(
            """
            Candidates are ranked by absolute mass error.
            The lowest ppm error is not automatically the
            chemically correct assignment.
            """,
            style={
                "fontWeight": "normal",
                "fontSize": "14px",
                "marginTop": "8px"
            }
        )
    ])

    # --------------------------------------------------------
    # RESULTS TABLE
    # --------------------------------------------------------

    table = dash_table.DataTable(

        data=df.to_dict(
            "records"
        ),

        columns=[
            {
                "name": column,
                "id": column
            }
            for column in df.columns
        ],

        sort_action="native",

        filter_action="native",

        page_action="native",

        page_size=20,

        style_table={
            "overflowX": "auto"
        },

        style_cell={
            "textAlign": "left",
            "padding": "8px",
            "fontFamily": "Arial",
            "fontSize": "14px"
        },

        style_header={
            "fontWeight": "bold",
            "backgroundColor": "#eeeeee"
        }
    )

    return message, table


# ============================================================
# RUN APP
# ============================================================

if __name__ == "__main__":

    app.run(
        debug=True
    )