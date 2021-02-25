import datetime
from typing import Tuple, List

import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import numpy as np
from pathlib import Path
import pandas as pd
from statsmodels.stats.proportion import proportion_confint

from base_classes import PlotConfig, GrowthNumbers, Dataset
import utils


INPUT_RAW_VIOLLIER = "../data/data_variants_per_day_and_grossregion_viollier.csv"
INPUT_RAW_RISCH = "../data/data_variants_per_day_risch.csv"
INPUT_CASES = "../data/data_cases_per_day_and_grossregion.csv"
OUTPUT_DIR = "../figures"

pconf = PlotConfig(
    colors=["#0D4A70", "#67916E", "#1883C6", "#99D9A4", "#73071F", "#5D8CF0"],
    greys=["#151624", "#8C8C8C", "#BFBFBF"]
)
global_generation_time = -1
r_viollier = 0.8298327047676343
r_risch = 0.8179411274559166
r_viollier_geneva = 0.8810484621340444


# --- Core functions ---

def calculate_growth_rates(
        data: pd.DataFrame,
        alpha: float,
        generation_time: float,
        reproduction_number: float
) -> GrowthNumbers:
    """
    Computes the growth rates from time series data assuming a logistic curve.

    :param data: columns: t:int, b117:int, original:int
    :param alpha:
    :param generation_time:
    :returns:
    """
    sequenced = data.b117 + data.original
    return utils.statsmodel_fit(alpha, data.t, data.b117, sequenced, generation_time, reproduction_number)


def draw_proportion_line_plot(
        ax: Axes,
        data: pd.DataFrame,
        generation_time,
        reproduction_number,
        start_date=datetime.datetime.strptime("2020-11-01", "%Y-%m-%d"),
        number_days=150,
        alpha=0.95,
        show_logistic_mle=True,
        show_logistic_ci=True,
        show_daily_points=True,
        show_daily_ci=True,
        show_legend=True
):
    """
    :param data: columns: date:datetime, year:int, week:int, t:int, region:str, cases:int, b117:int, original:int
    """
    growth = calculate_growth_rates(data, alpha, generation_time, reproduction_number)

    # Compute trajectories
    t_integer = np.array(range(number_days))
    t_date = [start_date + datetime.timedelta(days=x) for x in range(number_days)]
    model_pred = utils.statsmodel_predict(growth.statsmodel_model, t_integer, alpha)

    # Draw
    proportion = data.b117 / (data.b117 + data.original)
    daily_ci = proportion_confint(data.b117, data.b117 + data.original, alpha=1 - alpha, method="wilson")
    if show_logistic_mle:
        ax.plot(t_date, model_pred[0], label="Logistic fit", color=pconf.greys[0])
    if show_logistic_ci:
        ax.fill_between(t_date, model_pred[1], model_pred[0], color=pconf.greys[2], label="Confidence interval")
        ax.fill_between(t_date, model_pred[2], model_pred[0], color=pconf.greys[2])
    if show_daily_ci:
        yerr = np.vstack((proportion - daily_ci[0],
                          daily_ci[1] - proportion))
        ax.errorbar(data.date, proportion, yerr=yerr, fmt="none", c=pconf.colors[0])
    if show_daily_points:
        ax.scatter(data.date, proportion, s=14, c=pconf.colors[1],
                   label="Estimated daily proportion of B.1.1.7", zorder=10)
    if show_legend:
        ax.legend()
    ax.set_xlim(min(t_date), max(t_date))
    ax.set_ylim(0, 1)
    locator = mdates.AutoDateLocator(minticks=5, maxticks=10)
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)


def draw_model_vs_reality_plot(
        ax: Axes,
        data_confirmed_cases: pd.DataFrame,
        data_variants: pd.DataFrame,
        generation_time,
        model_original_r0,
        model_initial_cases=3200,
        model_initial_proportion=0.02,
        model_fitness_advantage=0.5,
        model_start_date=datetime.datetime.strptime("2021-01-01", "%Y-%m-%d"),
        model_end_date=datetime.datetime.strptime("2021-03-01", "%Y-%m-%d"),
        real_start_date=datetime.datetime.strptime("2021-01-01", "%Y-%m-%d"),
        real_sliding_window_size=7,
        show_legend=True,
        show_model_areas=True,
        show_reality_line_b117=True,
        show_reality_line_total=True,
        reality_line_b117_color=pconf.colors[0],
        reality_line_b117_label_suffix="",
        show_long_title=True,
        y_axis_lim=None
):
    """
    :param data_confirmed_cases: The number of confirmed cases per day, columns: date:datetime, cases:int
    :param data_variants: The number of B.1.1.7 vs. original, columns: date:datetime, b117:int, original:int
    """

    # Calculate R value
    model_b117_r0 = model_original_r0 * (1 + model_fitness_advantage)

    # Prepare time values
    number_days = (model_end_date - model_start_date).days + 1
    data_model = pd.DataFrame({
        "t": range(number_days)
    })
    data_model["date"] = pd.date_range(model_start_date, periods=number_days, freq="D")
    data_model["generation"] = data_model.t / generation_time

    # Calculate model's case numbers
    model_initial_variant = model_initial_cases * model_initial_proportion
    model_initial_original = model_initial_cases - model_initial_variant
    data_model["n_b117"] = model_initial_variant * (model_b117_r0 ** data_model.generation)
    data_model["n_original"] = model_initial_original * (model_original_r0 ** data_model.generation)

    # Join daily data with the real case numbers
    data_variants["proportion"] = data_variants.b117 / (data_variants.b117 + data_variants.original)
    daily_estimated_cases = data_confirmed_cases \
        .merge(data_variants, on="date", how="left") \
        .sort_values("date")
    daily_estimated_cases["b117"] = round(daily_estimated_cases.cases * daily_estimated_cases.proportion)
    daily_estimated_cases["original"] = round(daily_estimated_cases.cases - daily_estimated_cases.b117)
    daily_estimated_cases = daily_estimated_cases[["date", "b117", "original", "cases"]]
    daily_estimated_cases = daily_estimated_cases.rolling(real_sliding_window_size, center=True, on="date").mean()
    daily_estimated_cases = daily_estimated_cases[daily_estimated_cases.date >= real_start_date]

    # Draw plot
    if show_reality_line_total:
        ax.plot(daily_estimated_cases.date, daily_estimated_cases.cases,
                label="Data: Total number of confirmed cases (7-day average)", color=pconf.colors[1], linewidth=4)
    if show_reality_line_b117:
        ax.plot(daily_estimated_cases.date, daily_estimated_cases.b117,
                label="Data: Estimated number of B.1.1.7 cases (7-day average)" + reality_line_b117_label_suffix,
                color=reality_line_b117_color,
                linewidth=4)
    if show_model_areas:
        ax.stackplot(data_model.date, data_model.n_b117, data_model.n_original,
                     labels=["Model: Number of B.1.1.7 cases", "Model: Number of total cases"],
                     colors=[pconf.colors[2], pconf.colors[3]])
    if show_legend:
        ax.legend()
    if show_long_title:
        title = "Initial cases: {}, proportion of B.1.1.7: {}%\n" \
                "R-value of known variants: {}, fitness advantage: {}".format(
                    round(model_initial_cases), round(model_initial_proportion * 100, 1),
                    round(model_original_r0, 2), round(model_fitness_advantage, 2)
                )
        ax.set_title(title)
    if y_axis_lim is not None:
        ax.set_ylim(0, y_axis_lim)
    locator = mdates.AutoDateLocator(minticks=5, maxticks=10)
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)


# --- Data loading and transformation ---

def load_data_cases(file) -> pd.DataFrame:
    """
    :returns: columns: date:datetime, year:int, week:int, t:int, region:str, cases:int
    """
    data_cases = pd.read_csv(file)
    data_cases["date"] = pd.to_datetime(data_cases.date)
    data_cases["t"] = (data_cases.date - data_cases.date.min()).dt.days
    return data_cases


def load_data_variants(file) -> pd.DataFrame:
    """
    :returns: columns: date:datetime, region:str, sequenced:int, b117:int, original:int
    """
    data_variants = pd.read_csv(file)
    data_variants["date"] = pd.to_datetime(data_variants.date)
    data_variants["original"] = data_variants.sequenced - data_variants.b117
    return data_variants


def load_data(file_cases, file_variants) -> Dataset:
    return Dataset(
        load_data_cases(file_cases),
        load_data_variants(file_variants)
    )


def to_data_whole(data: Dataset) -> Dataset:
    """
    Removes region dimension by summing up the values for all regions on the same day
    """
    data_cases = data.cases \
        .groupby(["date", "year", "week", "t"]) \
        .aggregate({"cases": "sum"}) \
        .reset_index()
    data_variants = data.variants \
        .groupby(["date"]) \
        .aggregate({"sequenced": "sum", "b117": "sum", "original": "sum"}) \
        .reset_index()
    return Dataset(data_cases, data_variants)


def to_data_regions(data: Dataset) -> List[Tuple[str, Dataset]]:
    """
    :param data:
    :returns: a tuple: (region name, data)
    """
    result: List[Tuple[str, Dataset]] = []
    regions = pd.unique(data.cases.region)
    for region in regions:
        data_cases = data.cases[data.cases.region == region].reset_index(drop=True)
        data_variants = data.variants[data.variants.region == region].reset_index(drop=True)
        result.append((region, Dataset(data_cases, data_variants)))
    return result


def join_dataset(data: Dataset) -> pd.DataFrame:
    if "region" in data.cases:
        join_cols = ["date", "region"]
    else:
        join_cols = "date"
    return data.cases.merge(data.variants, on=join_cols, how="inner")


# -------------- Very paper-specific code --------------

def generate_fig_proportion_switzerland_viollier_risch() -> Figure:
    def _subplot(ax, data, start_date, reproduction_number):
        d = join_dataset(to_data_whole(data))
        draw_proportion_line_plot(ax, d, reproduction_number=reproduction_number,
                                  generation_time=global_generation_time, show_logistic_ci=True,
                                  show_daily_points=True, show_daily_ci=True,
                                  start_date=start_date, number_days=120)

    data_raw_viollier = load_data(INPUT_CASES, INPUT_RAW_VIOLLIER)
    data_raw_risch = load_data(INPUT_CASES, INPUT_RAW_RISCH)

    fig = plt.figure(figsize=(16, 6))

    ax = fig.add_subplot(1, 2, 1)
    _subplot(ax, data_raw_viollier, datetime.datetime.strptime("2020-12-14", "%Y-%m-%d"), r_viollier)
    ax.set_title("Switzerland - data from Viollier")

    ax = fig.add_subplot(1, 2, 2)
    _subplot(ax, data_raw_risch, datetime.datetime.strptime("2020-12-14", "%Y-%m-%d"), r_risch)
    ax.set_title("Switzerland - data from Risch")

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR + "/g" + str(global_generation_time) + "/fig_proportion_switzerland_viollier_risch.pdf",
                format="pdf")
    fig.savefig(OUTPUT_DIR + "/g" + str(global_generation_time) + "/fig_proportion_switzerland_viollier_risch.svg",
                format="svg")
    return fig


def generate_fig_proportion_regions() -> Figure:
    def _subplot(ax, data, start_date, show_legend):
        draw_proportion_line_plot(ax, data, reproduction_number=r_viollier,
                                  generation_time=global_generation_time, show_logistic_ci=True,
                                  show_daily_points=True, show_daily_ci=True,
                                  start_date=start_date, number_days=120, show_legend=show_legend)

    data_raw = load_data(INPUT_CASES, INPUT_RAW_VIOLLIER)

    fig = plt.figure(figsize=(16, 12))

    for i, (region, d) in enumerate(to_data_regions(data_raw)):
        ax = fig.add_subplot(4, 2, i + 1)
        _subplot(ax, join_dataset(d), start_date=datetime.datetime.strptime("2020-12-14", "%Y-%m-%d"),
                 show_legend=i == 0)
        ax.set_title(region)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR + "/g" + str(global_generation_time) + "/fig_proportion_regions.pdf", format="pdf")
    fig.savefig(OUTPUT_DIR + "/g" + str(global_generation_time) + "/fig_proportion_regions.svg", format="svg")
    return fig


def generate_fig_model_vs_reality_switzerland_viollier_risch(show_model=True) -> Figure:
    def _subplot(ax, d, reproduction_number, show_legend=False):
        growth = calculate_growth_rates(join_dataset(d), 0.95, global_generation_time, reproduction_number)
        data_confirmed_cases = d.cases[["date", "cases"]]
        data_variants = join_dataset(d)[["date", "t", "b117", "original"]]
        initial_cases = data_confirmed_cases[
            (data_confirmed_cases.date >= datetime.datetime.strptime("2020-12-29", "%Y-%m-%d")) &
            (data_confirmed_cases.date <= datetime.datetime.strptime("2021-01-04", "%Y-%m-%d"))]
        initial_cases = sum(initial_cases.cases) / 7
        draw_model_vs_reality_plot(ax, data_confirmed_cases, data_variants, model_original_r0=reproduction_number,
                                   generation_time=global_generation_time, model_initial_cases=initial_cases,
                                   model_fitness_advantage=growth.fc_mle, show_long_title=False,
                                   show_legend=show_legend, y_axis_lim=4000, show_model_areas=show_model)

    data_viollier = load_data(INPUT_CASES, INPUT_RAW_VIOLLIER)
    data_risch = load_data(INPUT_CASES, INPUT_RAW_RISCH)

    fig = plt.figure(figsize=(16, 6))

    ax = fig.add_subplot(1, 2, 1)
    _subplot(ax, to_data_whole(data_viollier), r_viollier, show_legend=True)
    ax.set_title("Switzerland - data from Viollier")

    ax = fig.add_subplot(1, 2, 2)
    _subplot(ax, to_data_whole(data_risch), r_risch)
    ax.set_title("Switzerland - data from Risch")

    fig.tight_layout()
    if show_model:
        fig.savefig(OUTPUT_DIR + "/g" + str(global_generation_time) + "/fig_model_vs_reality_switzerland_viollier_risch.pdf",
                    format="pdf")
    else:
        fig.savefig(OUTPUT_DIR + "/g" + str(global_generation_time) + "/fig_cases_switzerland_viollier_risch.svg",
                    format="svg")
    return fig


def generate_fig_cases_switzerland_joint() -> Figure:
    def _subplot(ax, d, reproduction_number, first=False):
        growth = calculate_growth_rates(join_dataset(d), 0.95, global_generation_time, reproduction_number)
        data_confirmed_cases = d.cases[["date", "cases"]]
        data_variants = join_dataset(d)[["date", "t", "b117", "original"]]
        initial_cases = data_confirmed_cases[
            (data_confirmed_cases.date >= datetime.datetime.strptime("2020-12-29", "%Y-%m-%d")) &
            (data_confirmed_cases.date <= datetime.datetime.strptime("2021-01-04", "%Y-%m-%d"))]
        initial_cases = sum(initial_cases.cases) / 7
        if first:
            draw_model_vs_reality_plot(ax, data_confirmed_cases, data_variants, model_original_r0=reproduction_number,
                                   generation_time=global_generation_time, model_initial_cases=initial_cases,
                                   model_fitness_advantage=growth.fc_mle, show_long_title=False,
                                   show_legend=False, y_axis_lim=4000, show_model_areas=False,
                                   reality_line_b117_label_suffix=", Viollier")
        else:
            draw_model_vs_reality_plot(ax, data_confirmed_cases, data_variants, model_original_r0=reproduction_number,
                                       generation_time=global_generation_time, model_initial_cases=initial_cases,
                                       model_fitness_advantage=growth.fc_mle, show_long_title=False,
                                       show_legend=True, y_axis_lim=4000, show_model_areas=False, show_reality_line_total=False,
                                       reality_line_b117_color=pconf.colors[2],
                                       reality_line_b117_label_suffix=", Risch")

    data_viollier = load_data(INPUT_CASES, INPUT_RAW_VIOLLIER)
    data_risch = load_data(INPUT_CASES, INPUT_RAW_RISCH)

    fig = plt.figure(figsize=(10, 7))

    ax = fig.add_subplot(1, 1, 1)
    _subplot(ax, to_data_whole(data_viollier), r_viollier, first=True)
    _subplot(ax, to_data_whole(data_risch), r_risch)
    ax.set_title("Switzerland")

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR + "/g" + str(global_generation_time) + "/fig_cases_switzerland_joint.svg",
                format="svg")
    return fig


def generate_fig_model_vs_reality_regions(show_model=True) -> Figure:
    def _subplot(ax, d, growth, show_legend=False):
        data_confirmed_cases = d.cases[["date", "cases"]]
        data_variants = join_dataset(d)[["date", "t", "b117", "original"]]
        initial_cases = data_confirmed_cases[
            (data_confirmed_cases.date >= datetime.datetime.strptime("2020-12-29", "%Y-%m-%d")) &
            (data_confirmed_cases.date <= datetime.datetime.strptime("2021-01-04", "%Y-%m-%d"))]
        initial_cases = sum(initial_cases.cases) / 7
        draw_model_vs_reality_plot(ax, data_confirmed_cases, data_variants, model_original_r0=r_viollier,
                                   generation_time=global_generation_time, model_initial_cases=initial_cases,
                                   model_fitness_advantage=growth.fc_mle, show_long_title=False,
                                   show_legend=show_legend, show_model_areas=show_model)

    data_raw = load_data(INPUT_CASES, INPUT_RAW_VIOLLIER)
    growth = calculate_growth_rates(join_dataset(to_data_whole(data_raw)), 0.95, global_generation_time, r_viollier)

    fig = plt.figure(figsize=(16, 12))

    for i, (region, d) in enumerate(to_data_regions(data_raw)):
        ax = fig.add_subplot(4, 2, i + 1)
        _subplot(ax, d, growth, show_legend=i == 0)
        ax.set_title(region)

    fig.tight_layout()
    if show_model:
        fig.savefig(OUTPUT_DIR + "/g" + str(global_generation_time) + "/fig_model_vs_reality_regions.pdf", format="pdf")
    else:
        fig.savefig(OUTPUT_DIR + "/g" + str(global_generation_time) + "/fig_cases_regions.svg", format="svg")
    return fig


def generate_tbl_growth_rates() -> pd.DataFrame:
    def _row(d, region, reproduction_number):
        g = calculate_growth_rates(d, 0.95, global_generation_time, reproduction_number)
        return {
            "Region": region,
            "$a$": "{:.2f} [{:.2f}; {:.2f}]".format(g.a_mle, g.a_ci.lower, g.a_ci.upper),
            "$t_0$": "{:.0f} [{:.0f}; {:.0f}]".format(g.t0_mle, g.t0_ci.lower, g.t0_ci.upper),
            "$f_d$": "{:.2f} [{:.2f}; {:.2f}]".format(g.fd_mle, g.fd_ci.lower, g.fd_ci.upper),
            "$f_c$": "{:.2f} [{:.2f}; {:.2f}]".format(g.fc_mle, g.fc_ci.lower, g.fc_ci.upper)
        }

    data_viollier = load_data(INPUT_CASES, INPUT_RAW_VIOLLIER)
    data_risch = load_data(INPUT_CASES, INPUT_RAW_RISCH)

    table = []

    # Switzerland
    table.append(_row(join_dataset(to_data_whole(data_viollier)), "Switzerland (Viollier)", r_viollier))
    table.append(_row(join_dataset(to_data_whole(data_risch)), "Switzerland (Risch)", r_risch))

    # Regions
    for region, d in to_data_regions(data_viollier):
        table.append(_row(join_dataset(d), region, r_viollier))

    df = pd.DataFrame(table)
    df.to_csv(OUTPUT_DIR + "/g" + str(global_generation_time) + "/tbl_growth_rates.csv", index=False)
    return df


def export_estimated_case_numbers_for_re():
    def _export(data, path):
        d = join_dataset(to_data_whole(data))
        d["sequenced"] = d.b117 + d.original
        d["b117"] = pd.to_numeric(round((d.b117 / d.sequenced) * d.cases), downcast="integer")
        d["original"] = pd.to_numeric(round((d.original / d.sequenced) * d.cases), downcast="integer")
        d[["date", "b117", "original"]].to_csv(path, index=False)

    # Switzerland
    _export(load_data(INPUT_CASES, INPUT_RAW_VIOLLIER), "Re/data/estimated_case_numbers_viollier.csv")
    _export(load_data(INPUT_CASES, INPUT_RAW_RISCH), "Re/data/estimated_case_numbers_risch.csv")

    # Regions
    for region, d in to_data_regions(load_data(INPUT_CASES, INPUT_RAW_VIOLLIER)):
        _export(d, "Re/data/estimated_case_numbers_viollier_" + region + ".csv")


def generate_fig_model_vs_reality_geneva():
    def _subplot(ax, d, growth, show_legend=False):
        data_confirmed_cases = d.cases[["date", "cases"]]
        data_variants = join_dataset(d)[["date", "t", "b117", "original"]]
        initial_cases = data_confirmed_cases[
            (data_confirmed_cases.date >= datetime.datetime.strptime("2020-12-29", "%Y-%m-%d")) &
            (data_confirmed_cases.date <= datetime.datetime.strptime("2021-01-04", "%Y-%m-%d"))]
        initial_cases = sum(initial_cases.cases) / 7
        draw_model_vs_reality_plot(ax, data_confirmed_cases, data_variants, model_original_r0=r_viollier_geneva,
                                   generation_time=global_generation_time, model_initial_cases=initial_cases,
                                   model_fitness_advantage=growth.fc_mle, show_long_title=False,
                                   show_legend=show_legend)

    data_raw = load_data(INPUT_CASES, INPUT_RAW_VIOLLIER)

    fig = plt.figure(figsize=(8, 3))

    for region, d in to_data_regions(data_raw):
        if region != "Lake Geneva region":
            continue
        growth = calculate_growth_rates(join_dataset(d), 0.95, global_generation_time, r_viollier_geneva)
        ax = fig.add_subplot(1, 1, 1)
        _subplot(ax, d, growth)
        ax.set_title(region)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR + "/g" + str(global_generation_time) + "/fig_model_vs_reality_geneva.pdf", format="pdf")
    return fig


def main():
    global global_generation_time
    font = {'size': 13}
    mpl.rc('font', **font)
    for g in [4.8, 5.2, 7.6]:
        global_generation_time = g
        Path(OUTPUT_DIR + "/g" + str(global_generation_time)).mkdir(parents=True, exist_ok=True)
        generate_fig_proportion_switzerland_viollier_risch()
        generate_fig_proportion_regions()
        generate_fig_model_vs_reality_switzerland_viollier_risch()
        generate_fig_model_vs_reality_regions()
        generate_fig_model_vs_reality_switzerland_viollier_risch(show_model=False)
        generate_fig_model_vs_reality_regions(show_model=False)
        generate_fig_cases_switzerland_joint()
        generate_tbl_growth_rates()
        generate_fig_model_vs_reality_geneva()
        if g == 4.8:
            export_estimated_case_numbers_for_re()


if __name__ == "__main__":
    main()
