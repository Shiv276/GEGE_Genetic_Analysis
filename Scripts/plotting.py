import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import plotly.graph_objects as go
from config import PROCESSED_DATA, FIGURES, GWAS_STRUCTURE_DIR


"""
The general plots used to inform a lot of decisions and interrogate data.
Not all will be used directly in the project or to even inform any decision.
"""

def load_processed_data():
    """
    Loading the pre-processed dfs.
    """
    df_ne = pd.read_csv(PROCESSED_DATA / "historical_ne.csv")
    froh = pd.read_csv(PROCESSED_DATA / "individual_froh.csv")
    froh_meta = pd.read_csv(PROCESSED_DATA / "individual_froh_metadata.csv")
    family_froh = pd.read_csv(PROCESSED_DATA / "family_froh.csv")
    pca = pd.read_csv(GWAS_STRUCTURE_DIR / "pca_with_metadata.csv")

    return df_ne, froh, froh_meta, family_froh, pca



def plot_historical_ne(df_ne, save=True):
    """
    Plots the effective population size (ne) over
    generations as provided by Davidson., et al (2026)

    NOTE* The save flag only exists for diagnostic testing.
    """
    fig, ax = plt.subplots(figsize=(7, 4))

    ax.plot(df_ne["Generations_ago"], df_ne["Ne"], marker="o")

    ax.set_xlabel("Generations ago")
    ax.set_ylabel(r"Effective population size ($N_e$)")

    #Makes time increase on the x axis (right-most would be current day)
    ax.invert_xaxis()


    # Approximate timing used by the original authors for Bd arrival
    ax.axvline( 10, linestyle="--", alpha=0.7, label="~Bd arrival")
    ax.legend(loc="upper left", frameon=False)


    fig.tight_layout()
    if save:
        fig.savefig(FIGURES / "historical_effective_population_size.png", dpi=300, bbox_inches="tight")

    return fig, ax



#Individual F_ROH Distribution
def plot_froh_distribution(froh, save=True):

    fig, ax = plt.subplots(figsize=(7, 4))

    ax.hist(froh["F_ROH_pct"], bins=30)

    mean_froh = froh["F_ROH_pct"].mean()
    median_froh = froh["F_ROH_pct"].median()

    #Dummy lines to create label handles for a legend to report mean/median
    mean_label = ax.plot([], [], linestyle="none", label=f"Mean = {mean_froh:.2f}%")
    median_label = ax.plot([], [], linestyle="none", label=f"Median = {median_froh:.2f}%")
    ax.legend(frameon = False)

    ax.set_xlabel("Genome contained in runs of homozygosity (%)")
    ax.set_ylabel("Number of individuals")



    fig.tight_layout()

    if save:
        fig.savefig(FIGURES / "individual_froh_distribution.png", dpi=300, bbox_inches="tight")
    else:
        pass

    return fig, ax




#Family level ROH
def plot_family_froh(family_froh, save=True):
    """
    Plots the ROH distribution across the 54 families.
    """

    #Sorting so that the plot looks cleaner.
    family_froh = family_froh.sort_values("mean_FROH_pct").copy()

    fig, ax = plt.subplots(figsize=(7, 10))

    ax.errorbar(
        family_froh["mean_FROH_pct"],
        family_froh["clutch"].astype(str),
        xerr=family_froh["se_FROH_pct"],
        fmt="o",
        capsize=2
    )

    ax.set_xlabel(r"Mean $F_{ROH}$ (%) ± SE")
    ax.set_ylabel("Pedigree family")

    fig.tight_layout()

    if save:
        fig.savefig(FIGURES / "family_froh.png", dpi=300, bbox_inches="tight")

    return fig, ax



def plot_roh_number_vs_froh(froh, save=True):
    """
    Plots ROH count against F_ROH to assess whether variation in genomic homozygosity
    is primarily associated with differences in ROH number. 

    Mostly a small exploratory check that F_ROH is largely associated with ROH number
    """
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.scatter(froh["number_ROH"], froh["F_ROH_pct"], alpha=0.6)

    ax.set_xlabel("Number of runs of homozygosity")

    ax.set_ylabel(r"$F_{ROH}$ (%)")

    fig.tight_layout()

    if save:
        fig.savefig(FIGURES / "roh_number_vs_froh.png", dpi=300, bbox_inches="tight")

    return fig, ax


def plot_pca_by_zoo(pca, save=True):

    #Create a template for plotting
    fig, ax = plt.subplots(figsize=(6, 5))

    #Groupby zoo
    groups = pca.groupby("zoo")

    for zoo, group in groups:
        ax.scatter(group["PC1"], group["PC2"], label=zoo, alpha=0.7, s=20)

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")

    ax.legend(title="Zoo")
    ax.invert_yaxis()
    fig.tight_layout()

    if save:
        fig.savefig(
            FIGURES /
            "pca_pc1_pc2_by_zoo.png",
            dpi=300,
            bbox_inches="tight"
        )

    return fig, ax



def plot_pca_family_slider(pca):
    """
    Attempt at an interactive PCA plot that shows variance among each family/clutch. 

    Aiming to turn this into an animation / gif to visually describe varaince across the 54 populations step-wise.
    Method heavily informed / helped from AI (was very difficult to figure the loop structure and new plotly tools)
    """

    fig = go.Figure()

    #Prepare the families series for showing each family per scroll.
    families = sorted(pca["clutch"].unique())

    #Adding a plotly graph object (scatter plot) as a trace (visual data layers for GOs)
    #These represent the background frogs / points that aren't currently in focus (hence grey and opaque).
    fig.add_trace(go.Scatter( x=pca["PC1"], y=pca["PC2"], mode="markers", marker=dict(size=5, color="lightgrey", opacity=0.2)))
  
    #Looping through families in order to create a unique trace state for each so that their relative sparsity can be highlighted
    for i, family in enumerate(families):

        #Filters the df for only frogs that are in the current iterations family
        group = pca[pca["clutch"] == family]

        #Adds a new scatter (trace) to overlay the current one, representing the active state for this family
        fig.add_trace(go.Scatter(x=group["PC1"], y=group["PC2"], mode="markers", marker=dict(size=9, color='MediumPurple'), name=f"Family {family}", 
                                 visible = (i == 0) #Only the first trace is visible at the beginning
                                 ))
                                 

    #Create step controls so plotly can show the correct trace when slider is moved.
    steps = []

    #Loop through each family once more.
    for i, family in enumerate(families):

        #First value belongs to background trace. Background should always stay visible.
        visible = [True]

        #Loop through all 54 families again in an inner loop.
        #Visible[0] remains True
        for j in range(len(families)):

            visible.append(j == i)
            #THEN, the visible list gets built with a bunch of False's EXCEPT where the current trace is (i).
            #So if i == 4 (family 5), visible becomes: [True (background), False (fam1), False (fam2), False (fam3), False (fam4), TRUE (fam5), False (fam6), ..., False (fam54)]
            #Only the first (background) and FIFTH traces have a True value in its visible list.
            
            #The loop resets after and the visible list resets to [True] with it.
    
        step = dict(
            method="update", #Instructions to update the existing figure
            args=[{"visible": visible},  {"title": f"PCA — Pedigree family {family}"}], 
            #visible tells plotly which trace to show/hide based on the visible list just built
            #Title also gets updated based on the current family number so viewers can see the progression across groups.
            label=str(family) #Determines the text that appears on the slider (being the current family number)
        )

        steps.append(step) #Adds that method/function into the steps list so that the slider can perform them when selected.
        #So steps would accumulate all 54 instructions (dictionaries) by the end of the outer loop.



    #Add slider and format figure

    fig.update_layout(
        title=f"PCA — Pedigree family {families[0]}",
        xaxis_title="PC1",
        yaxis_title="PC2",
        template="plotly_white",
        showlegend=False,
        sliders=[dict(active=0, currentvalue=dict(prefix="Family: "), steps=steps)], #Sliders is called on a list because we can have multiple. One is used here. (list of dicts)
        
        #active is the initially selected slider step (starts at 0)
        #currentvalue controls the text displayed alongside the slider indicating its current selection. Its combined with the label of the selected step (line 227)
        #steps as a parameter tells the slider what to do. We set it to the .update() method in our steps list (somewhat poorly named in retrospect, but gets the idea across).
        width=900,
        height=650
    )

    """
    As a step is updated, plotly effectively does (assuming active=3):

    Plotly.update(
    graph,
    {"visible": [true, false, false, true, false]},
    {"title": "PCA - Pedigree family 3"}
    );
    """

    fig.show()










def main():

    df_ne, froh, froh_meta, family_froh, pca = load_processed_data()

    #plot_historical_ne(df_ne)
    #plot_froh_distribution(froh)
    #plot_family_froh(family_froh)
    #plot_roh_number_vs_froh(froh)
    #plot_pca_by_zoo(pca, save=False)
    #plot_pca_family_centroids(pca, save=False)
    #plot_pca_by_family(pca, save=False)
    #plot_interactive_pca_family_slider(pca)
    plot_pca_family_slider(pca)
    plt.show()


if __name__ == "__main__":
    main()