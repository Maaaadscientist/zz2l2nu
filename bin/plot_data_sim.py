#!/usr/bin/env python

"""Plots distributions of data and simulation."""

import argparse
import math
import os

import numpy as np

import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import yaml

from pyroothist import Hist1D


def fill_histograms(samples, selections, variables):
    histograms = {}
    for isample, sample in enumerate(samples):
        chain = ROOT.TChain('Vars')
        for path in sample['files']:
            chain.AddFile(path)
        data_frame = ROOT.RDataFrame(chain)
        proxies = {}

        for iselection, selection in enumerate(selections):
            df_filtered = data_frame.Filter(selection['formula']).Define(
                '_weight', selection['weight'])
            for ivariable, variable in enumerate(variables):
                bin_def = variable['binning'][selection['tag']]
                r = bin_def['range']
                num_bins = bin_def['num_bins']

                if variable.get('scale', 'linear') == 'log':
                    binning = np.geomspace(r[0], r[1], num_bins + 1)
                else:
                    binning = np.linspace(r[0], r[1], num_bins + 1)
                
                proxies[iselection, ivariable] = df_filtered.Define(
                    '_var', variable['formula']
                ).Histo1D(
                    ROOT.RDF.TH1DModel('', '', num_bins, binning),
                    '_var', '_weight')

        for key, proxy in proxies.items():
            histograms[isample, key[0], key[1]] = Hist1D(proxy.GetValue())

    return histograms


def plot_data_sim(variable_info, data_hist, sim_hists_infos,
                  output='fig.pdf', selection_label=''):
    fig = plt.figure()
    fig.patch.set_alpha(0.)
    gs = mpl.gridspec.GridSpec(
        3, 2, hspace=0., wspace=0.,
        height_ratios=[2, 1, 1], width_ratios=[6, 1]
    )
    axes_distributions = fig.add_subplot(gs[0, 0])
    axes_composition = fig.add_subplot(gs[1, 0])
    axes_residuals = fig.add_subplot(gs[2, 0])

    sim_hist_total = Hist1D()
    for s in sim_hists_infos:
        sim_hist_total += s[0]

    binning = data_hist.binning
    widths = binning[1:] - binning[:-1]

    handle_sim = axes_distributions.hist(
        binning[:-1], weights=sim_hist_total.contents[1:-1] / widths,
        bins=binning, histtype='stepfilled',
        ec='dodgerblue', fc=mpl.colors.to_rgba('dodgerblue', 0.2)
    )[2][0]
    handle_data = axes_distributions.errorbar(
        (binning[:-1] + binning[1:]) / 2, data_hist.contents[1:-1] / widths,
        yerr=data_hist.errors[1:-1] / widths,
        marker='o', color='black', ls='none'
    )
    axes_distributions.set_yscale('log')

    axes_composition.hist(
        [binning[:-1]] * len(sim_hists_infos), bins=binning,
        weights=[
            s[0].contents[1:-1] / sim_hist_total.contents[1:-1]
            for s in reversed(sim_hists_infos)
        ],
        color=[s[1]['color'] for s in reversed(sim_hists_infos)],
        stacked=True
    )
    axes_composition.set_ylim(0., 1.)
    axes_composition.yaxis.set_major_locator(mpl.ticker.NullLocator())
    axes_composition.yaxis.set_minor_locator(mpl.ticker.NullLocator())

    sim_errors = np.empty(len(binning))
    sim_errors[:-1] = sim_hist_total.errors[1:-1] \
        / sim_hist_total.contents[1:-1]
    sim_errors[-1] = sim_errors[-2]  # Needed for fill_between
    axes_residuals.fill_between(
        binning, 1. + sim_errors, 1. - sim_errors,
        step='post', lw=0., color='0.75'
    )
    axes_residuals.errorbar(
        (binning[:-1] + binning[1:]) / 2,
        data_hist.contents[1:-1] / sim_hist_total.contents[1:-1],
        yerr=data_hist.errors[1:-1] / sim_hist_total.contents[1:-1],
        marker='o', color='black', ls='none'
    )
    axes_residuals.axhline(1., ls='dashed', lw=0.8, color='black')
    axes_residuals.set_ylim(0., 2.)

    for axes in [axes_distributions, axes_composition, axes_residuals]:
        axes.set_xlim(binning[0], binning[-1])

    if variable_info.get('scale', 'linear') == 'log':
        for axes in [axes_distributions, axes_composition, axes_residuals]:
            axes.set_xscale('log')

        # Provide a formatter for minor ticks so that they get labelled.
        # Also change the formatter for major ticks in order to get a
        # consistent formatting (1000 intead of 10^3).
        axes_residuals.xaxis.set_major_formatter(mpl.ticker.LogFormatter())
        axes_residuals.xaxis.set_minor_formatter(
            mpl.ticker.LogFormatter(minor_thresholds=(2, 0.4)))

    for axes in [axes_distributions, axes_composition]:
        #axes.set_xticklabels([])
        axes.xaxis.set_major_formatter(mpl.ticker.NullFormatter())
        axes.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

    if 'unit' in variable_info:
        xlabel = '{} [{}]'.format(
            variable_info['label'], variable_info['unit'])
        ylabel_distributions = r'$\langle\mathrm{{Events}}\//\/' \
            r'\mathrm{{{}}}\rangle$'.format(variable_info['unit'])
    else:
        xlabel = variable_info['label']
        ylabel_distributions = r'$\langle\mathrm{Events}\//\/' \
            r'\mathrm{unit}\rangle$'

    axes_residuals.set_xlabel(xlabel)
    axes_distributions.set_ylabel(ylabel_distributions)
    axes_composition.set_ylabel('Compos.')
    axes_residuals.set_ylabel(r'$\mathrm{Data}\//\/\mathrm{Exp.}$')
    fig.align_ylabels()

    if selection_label:
        axes_distributions.text(
            1., 1.002, selection_label,
            ha='right', va='bottom', transform=axes_distributions.transAxes
        )

    # Manually construct a common legend for all panels
    handles=[handle_data, handle_sim]
    labels=['Data', 'Exp.']

    for s in sim_hists_infos:
        info = s[1]
        handles.append(mpl.patches.Patch(color=info['color']))
        labels.append(info['label'])

    axes_distributions.legend(
        handles, labels,
        bbox_to_anchor=(1., 1.), loc='upper left', frameon=False
    )

    fig.savefig(output)
    plt.close(fig)


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser(description=__doc__)
    arg_parser.add_argument('config', help='Configuration file')
    arg_parser.add_argument(
        '-o', '--output', default='fig',
        help='Directory for produced figures')
    arg_parser.add_argument(
        '-y', '--year', default='2016',
        help='Year of data taking')
    args = arg_parser.parse_args()


    with open(args.config) as f:
        config = yaml.safe_load(f)

    selections = config['selections']
    variables = config['variables']
    
    data_sample = config['samples']['data']
    sim_samples = config['samples']['simulation']

    file_prefix = config['samples'].get('path_prefix', '')
    for sample in [data_sample] + sim_samples:
        sample['files'] = [file_prefix + p for p in sample['files']]

    for selection in selections:
        try:
            os.makedirs(os.path.join(args.output, selection['tag']))
        except FileExistsError:
            pass

    # plt.style.use('my-standard')

    
    histograms = fill_histograms(
        [data_sample] + sim_samples, selections, variables)

    # Include under- and overflows into neighbouring bins and clip
    # negative bin contents
    for hist in histograms.values():
        hist.contents[1] += hist.contents[0]
        hist.contents[-2] += hist.contents[-1]
        hist.contents[0] = hist.contents[-1] = 0.

        hist.errors[1] = math.hypot(hist.errors[0], hist.errors[1])
        hist.errors[-2] = math.hypot(hist.errors[-2], hist.errors[-1])
        hist.errors[0] = hist.errors[-1] = 0.

        np.clip(hist.contents, 0., None, out=hist.contents)


    for ivariable, variable in enumerate(variables):
        for iselection, selection in enumerate(selections):
            plot_data_sim(
                variable,
                histograms[0, iselection, ivariable],
                [
                    (histograms[i + 1, iselection, ivariable], sim_samples[i])
                    for i in range(len(sim_samples))
                ],
                output=os.path.join(args.output, selection['tag'],
                                    variable['tag'] + '.pdf'),
                selection_label='{}, {}'.format(selection['label'], args.year)
            )

