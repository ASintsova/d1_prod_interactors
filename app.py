from time import sleep
from dash import Dash, dcc, html, dash_table, Input, Output, State, ctx
import dash_bootstrap_components as dbc
import pandas as pd
import webbrowser
from dash.exceptions import PreventUpdate
import numpy as np
import plotly.graph_objects as go
from dash_bootstrap_templates import load_figure_template
import requests
from dash.dash_table.Format import Format, Scheme


app = Dash(
    __name__,
    external_stylesheets=[dbc.themes.SPACELAB, dbc.icons.FONT_AWESOME],
)
server = app.server

load_figure_template("spacelab")

# Load the data
lfc_data = pd.read_csv("data/proteomics_lfc.csv.gz")
lfc_data['UniProt'] = (lfc_data["Uniprot.ID"]
                       .str.split("-", expand=True)[0]
                       .str.strip())
lfc_data['link'] = (lfc_data['UniProt']
                    .apply(lambda x: f'https://www.uniprot.org/uniprotkb/{x}'))

lfc_data['UniProt'] = lfc_data['UniProt'].apply(lambda x: f"[{x}](https://www.uniprot.org/uniprotkb/{x})")
lfc_data['-log10 (padj)'] = -1 * np.log10(lfc_data['pval'])
lfc_data = lfc_data.rename(columns={'Gene name': 'Gene'})
to_show = ['UniProt', 'Gene', 'Tissue', 'Bait', 'LFC', 'pval']
embryo_data = lfc_data[lfc_data.Tissue == "Embryo"]
ovary_data = lfc_data[lfc_data.Tissue == "Ovary"]
testis_data = lfc_data[lfc_data.Tissue == "Testis"]


# App Layout

# Data table
lfc_data_table = dash_table.DataTable(
    id='lfc_data_table',
    # columns=[
    #     {"name": i, "id": i, "deletable": True, "selectable": True} for i in to_show
    # ],
    columns=[
        dict(id='UniProt', name='UniProt',presentation= 'markdown'),
        dict(id='Gene', name='Gene'),
        dict(id='Tissue', name='Tissue'),
        dict(id='Bait', name='Bait'),
        dict(id='LFC', name='LFC', type='numeric',
             format=Format(precision=2, scheme=Scheme.fixed)),
        dict(id='pval', name='pval', type='numeric',
             format=Format(precision=2, scheme=Scheme.exponent))

    ],
    data=lfc_data.to_dict('records'),
    editable=False,
    filter_action="native",
    sort_action="native",
    sort_mode="multi",
    column_selectable="single",
    selected_columns=[],
    selected_rows=[],
    page_action="native",
    page_current=0,
    page_size=10,
    style_cell={
                'overflow': 'hidden',
                'textOverflow': 'ellipsis',
                'minWidth': '100px', 'width': '120px', 'maxWidth': '120px',
                'padding': '5px',
                
    },
    style_table={'height': '30vw', 'overflowY': 'auto'},
    style_header={
        'backgroundColor': 'white',
        'fontWeight': 'bold',
        'fontsize': 8,
        'font-family': 'sans-serif'
    },
    style_data={'fontsize': 6, 'font-family': 'sans-serif'},
    
)


intro_card = html.Div(

    [
        html.H5("Welcome", className="card-title my-1"),
        dcc.Markdown("Here you can browse through the data from our most recent study (Chavan, Skrutl et al., 2023). \
                         We performed quantitative MS using 3 baits (D1, Prod and Piwi) to identify putative interactors \
                         (L2FC>1, p<0.05) across Drosophila embryos, ovaries and testes. See the **About** section for more info on how to use this app.",
                     className="mx-0 my-0"),
        dbc.Container([
        dbc.Button("Go to the preprint", color="primary", href="https://www.biorxiv.org/content/10.1101/2023.07.11.548599v1", target="_b",
                   class_name="p-2 my-0"), 
        dbc.Button(html.Span([html.I(className="fab fa-github mx-1"), "View on GitHub"]), color="secondary", href="https://github.com/ASintsova/d1_prod_interactors", target="_b", 
                   className="p-2 pr-3 mx-4 my-1 ps-2 pe-3", )]),
    ], className="m-2"

)


footnote_card = html.Div(

    [
        html.H5("About", className="card-title"),
        dcc.Markdown("""On the left side is a table that allows you to filter hits based on gene, tissue, bait, LFC and p-value. 
                         For example, you can type `>1` in the **LFC** column to filter hits that have a log2FC > 1. You can also type `D1` in the **Bait** column 
                         to select for D1-specific interactors, while typing `ne Prod` will selectively exclude Prod-specific interactions. The UniProt column links each
                         protein to its respective UniProt page. The filtered results can be submitted to [STRING-db](https://string-db.org/) for functional analyses. 
                         The right side shows the results as volcano plots across different tissues. Clicking on individual point will show information about this protein
                         in the table on the left. All the original pulldown data can be found in Tables S1-3 and Table S6 from the preprint. 
                         Information about DNA repair and transposon repression proteins that were pulled down by D1 or Prod across the different samples can be found in Table S4 and S5 respectively.
"""),

    ], className="m-1 mx-0"

)

# Tissue specific


def get_bait_card(tissue_df, tissue):
    bait_selection_card = dbc.RadioItems(
        id=f"{tissue}_bait",
        options=tissue_df.Bait.unique(),
        value=tissue_df.Bait.unique()[0],
        inline=True,
        className="p-1 m-2 border rounded ",
    )
    return bait_selection_card


def get_gene_selection_card(tissue_df, tissue):
    gene_selection_card = dcc.Dropdown(
        id=f'{tissue}_gois',
        options=tissue_df['Gene'].unique(),
        multi=True,
        className="m-0 p-2"
    )
    return gene_selection_card


# Embryo
volc_graph_embryo = dcc.Graph(
    id='volc_embryo', className="m-0"), html.Pre(id='embryo_data')
embryo_bait_selection_card = get_bait_card(embryo_data, 'Embryo')
embryo_gene_selection_card = get_gene_selection_card(embryo_data, 'Embryo')

# Ovary
volc_graph_ovary = dcc.Graph(
    id='volc_ovary', className="m-0"), html.Pre(id='ovary_data')
ovary_bait_selection_card = get_bait_card(ovary_data, 'Ovary')
ovary_gene_selection_card = get_gene_selection_card(ovary_data, 'Ovary')

# Testis
volc_graph_testis = dcc.Graph(
    id='volc_testis', className="m-0"), html.Pre(id='testis_data')
testis_bait_selection_card = get_bait_card(testis_data, 'Testis')
testis_gene_selection_card = get_gene_selection_card(testis_data, 'Testis')

tab1 = dbc.Tab([dbc.Row([dbc.Col(embryo_bait_selection_card, lg=6,
                                 className="pt-2"),
                         dbc.Col(embryo_gene_selection_card, lg=6,
                                 className="pt-2")], className='m-0 p-0'),
                dbc.Row(volc_graph_embryo, className='m-0 h-100'),
                ],
               tab_id="tab-1",
               label="Embryo",

               )

tab2 = dbc.Tab([dbc.Row([dbc.Col(ovary_bait_selection_card, lg=6,
                                 className="pt-2"),
                         dbc.Col(ovary_gene_selection_card, lg=6,
                                 className="pt-2")], className='m-0 p-0'),
                dbc.Row(volc_graph_ovary, className='m-0 h-100'),
                ],
               tab_id="tab-2",
               label="Ovary",

               )
tab3 = dbc.Tab([dbc.Row([dbc.Col(testis_bait_selection_card, lg=6,
                                 className="pt-2"),
                         dbc.Col(testis_gene_selection_card, lg=6,
                                 className="pt-2")], className='m-0 p-0'),
                dbc.Row(volc_graph_testis, className='m-0 h-100'),
                ],
               tab_id="tab-3",
               label="Testis",

               )


tabs = dbc.Tabs(
    [tab1,
     tab2,
        tab3,
     ],
    id="tabs",
    active_tab="tab-1",
    className="mt-2",
)

app.layout = dbc.Container([
    dbc.Row(
        dbc.Col(
            html.H2(
                "Multi-tissue proteomics of satellite DNA-binding proteins",
                className="text-center bg-primary text-white p-2",
            ),
        )

    ),
    dbc.Row(
        [
            dbc.Col([intro_card, lfc_data_table,
                    html.Div(id="test_mrk"),
                    dbc.Container(
                     [dbc.Button('Generate PPI network for selected genes via STRING-db',
                                id='submit-val', n_clicks=0, color="primary", className="ms-1 me-3"),
                                dbc.Button('Clear filters', id='clear-filter', color="secondary",
                                           className="me-1")]),
                     html.Div(id='string_link', className="ms-3")], width=12, lg=6, className="m-0 p-0"),

            dbc.Col([tabs],

                    width=12,
                    lg=6,
                    className="pt-4",
                    ),


        ],
        className="ms-1",
    ),
    dbc.Row(
        dbc.Col(
            [footnote_card], className="m-3"
        )),

], fluid=True, )


def proteomics_volcano(df, gois):
    df['hits'] = (abs(df.LFC) > 1) & (df.pval < 0.05)
    df['colors'] = ['#446e9b' if i else '#999' for i in df.hits.values]
    text = [f'Gene: {string1}<br>Uniprot: {string2}<br>LFC: {string3}<br>pval: {string4}'
            for string1, string2, string3, string4 in zip(df['Gene'], df['Uniprot.ID'], df['LFC'], round(df['pval'], 7))]
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df.LFC, y=df["-log10 (padj)"],
                             hoverinfo='text',
                             text=text,
                             customdata=df.UniProt,
                             mode='markers',
                             marker=dict(color=df.colors)))

    fig.add_hline(y=1.3, fillcolor='grey', line_width=0.5, line_dash="dash")
    fig.add_vline(x=-1, fillcolor='grey', line_width=0.5, line_dash="dash")
    fig.add_vline(x=1, fillcolor='grey', line_width=0.5, line_dash="dash")

    if gois:
        gois_df = df[df['Gene'].isin(gois)]
        if not gois_df.empty:
            text = [f'Gene: {string1}<br>Uniprot: {string2}<br>LFC: {string3}<br>pval: {string4}'
                    for string1, string2, string3, string4 in zip(gois_df['Gene'],
                                                         gois_df['Uniprot.ID'], gois_df['LFC'], gois_df['pval'])]
            fig.add_trace(go.Scatter(x=gois_df.LFC, y=gois_df["-log10 (padj)"],
                                     hoverinfo='text',
                                     text=text,
                                     customdata=gois_df.UniProt,
                                     mode='markers',
                                     marker=dict(
                                         color='#d47500',
                                         size=20,
                                         line=dict(
                                             color='black',
                                             width=1
                                         )
            ),
            ))

    fig.update_layout(showlegend=False,  xaxis_title='LFC', yaxis_title='-log10 (pval)',
                      height=500,
                      margin={
                                        'l': 10, 'b': 0, 't': 20, 'r': 5
                                    })
    return fig


@app.callback(
  [Output("lfc_data_table", "data"),
   Output("lfc_data_table", "filter_query")],
  [Input("clear-filter", "n_clicks"),
   Input("volc_embryo", "clickData"),
   Input("volc_ovary", "clickData"),
   Input("volc_testis", "clickData")],

)
def resetFilter(n_clicks, clickData, clickData_ov, clickData_test):
    triggered = ctx.triggered_id
    if triggered == 'clear-filter':
        return lfc_data.to_dict('records'), ""
    elif triggered == 'volc_embryo': 
        filtered = lfc_data[lfc_data['UniProt'] == clickData['points'][0]['customdata']]
        return filtered.to_dict('records'), ""
    elif triggered == 'volc_ovary':
        filtered = lfc_data[lfc_data['UniProt'] == clickData_ov['points'][0]['customdata']]
        return filtered.to_dict('records'), ""
    elif triggered == 'volc_testis':
        filtered = lfc_data[lfc_data['UniProt'] == clickData_test['points'][0]['customdata']]
        return filtered.to_dict('records'), ""
    else:
        raise PreventUpdate
      



# @app.callback(
#     Output("lfc_data_table", "data"),
#     Input("volc_embryo", "clickData"),
#     prevent_initial_call=True,
# )
# def verify(clickData):
#     if clickData:
#         print("Verify triggered")
#         # those are the rownumbers in users_to_verify 
#         # if users_to_verify is not dynamic, just filter
#         # it by row number and grab user_ids or other cols
#         filtered = lfc_data[lfc_data['UniProt'] == clickData['points'][0]['customdata']]
#         print(filtered.index)
#         print(filtered.to_dict('records'))
#         return filtered.to_dict('records')

#     else:
#         return lfc_data.to_dict('records')

    


# @app.callback(
#     [Output("lfc_data_table", "data")],
#     [Input("volc_embryo", "clickData")]
# )
# def on_trace_click(clickData):
#     """Listen to click events and update table, passing filtered rows"""
#     if clickData:
#         p = clickData["points"][0]

#         # here, use 'customdata' property of clicked point, 
#         # could also use 'curveNumber', 'pointIndex', etc.
#         if 'customdata' in p:
#             key = p['customdata']
#         df_f = get_corresponding_rows(lfc_data, key)
#         return df_f.to_dict('records')
#     else:
#         raise PreventUpdate

# def get_corresponding_rows(df, my_key):
#     """Filter df, return rows that match my_key"""
#     return df[df['UniProt'] == my_key]
#return filtered_df.sort_values(by=['lastUpdated']).to_dict('records'), [row_id]



# @app.callback(
#     Output('embryo_data', 'children'),
#     Input('volc_embryo', 'clickData'))
# def open_embryo_url(clickData):
#     if clickData:
#         webbrowser.open(clickData["points"][0]["customdata"])
#     else:
#         raise PreventUpdate


@app.callback(
    Output("volc_embryo", "figure"),
    Input("Embryo_bait", "value"),
    Input('Embryo_gois', 'value')
)
def update_volcano_embryo(embryo_bait, gois):
    df = lfc_data[(lfc_data.Tissue == 'Embryo') & (
        lfc_data.Bait == embryo_bait)].copy()
    fig = proteomics_volcano(df, gois)
    return fig


@app.callback(
    Output('ovary_data', 'children'),
    Input('volc_ovary', 'clickData'))
def open_ovary_url(clickData):
    if clickData:
        webbrowser.open(clickData["points"][0]["customdata"])
    else:
        raise PreventUpdate


@app.callback(
    Output("volc_ovary", "figure"),
    Input("Ovary_bait", "value"),
    Input('Ovary_gois', 'value')
)
def update_volcano_ovary(embryo_bait, gois):
    df = lfc_data[(lfc_data.Tissue == 'Ovary') & (
        lfc_data.Bait == embryo_bait)].copy()
    fig = proteomics_volcano(df, gois)
    return fig


@app.callback(
    Output('testis_data', 'children'),
    Input('volc_testis', 'clickData'))
def open_testis_url(clickData):
    if clickData:
        webbrowser.open(clickData["points"][0]["customdata"])
    else:
        raise PreventUpdate


@app.callback(
    Output("volc_testis", "figure"),
    Input("Testis_bait", "value"),
    Input('Testis_gois', 'value')
)
def update_volcano_testis(embryo_bait, gois):
    df = lfc_data[(lfc_data.Tissue == 'Testis') & (
        lfc_data.Bait == embryo_bait)].copy()
    fig = proteomics_volcano(df, gois)
    return fig


@app.callback(
    Output('test_mrk', 'children'),
    Input('lfc_data_table', "derived_virtual_data"),
    # prevent_initial_call=True
)
def link_to_string(rows):
    df = pd.DataFrame(rows)
    if not df.empty:
        gene_names = df['Gene'].unique()
        uniprot_names = df['Uniprot.ID'].unique()
        return html.Div([dcc.Markdown(f"Number of unique gene names: {len(gene_names)}. Number of unique UNIPROT IDs: {len(uniprot_names)}"
                                      ), ], className="m-2")
    return html.Div([dcc.Markdown("No genes selected")])


@app.callback(
    Output('string_link', 'children'),
    Input('submit-val', 'n_clicks'),
    State('lfc_data_table', "derived_virtual_data"),
    prevent_initial_call=True
)
def link_to_string2(n_clicks, rows):

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = 'tsv-no-header'
    method = 'get_link'

    request_url = "/".join([string_api_url, output_format, method])
    species = 7227
    df = pd.DataFrame(rows)
    gene_names = [gene for gene in df['Gene'].unique() if gene]
    if len(gene_names) < 550:
        params = {
            "identifiers": "\r".join(gene_names),  # your protein
            "species": species,  # species NCBI identifier
            "network_flavor": "confidence",  # show confidence links
            "caller_identity": "explodata"  # your app name
        }
        network = requests.post(request_url, data=params)
        network_url = network.text.strip()
        link = html.Div([dbc.Button(
            f"View PPI network for {len(gene_names)} genes", href=network_url, 
            target="_blank", color="success", class_name="p-2 my-2")])
        sleep(0.5)
    else:
        link = html.Div([dbc.Alert("Too many genes selected", class_name="p-2 my-2 w-50", color="warning")])
    # if st_col.button('Get STRING network'):
    return link


if __name__ == "__main__":
    app.run_server(debug=True)
