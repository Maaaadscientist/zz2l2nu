#!/usr/bin/env python

"""Plots systematic variations from file with templates."""

import argparse
import os

import matplotlib as mpl
from matplotlib import pyplot as plt

import ROOT

from hzz import Hist1D, mpl_style

ROOT.PyConfig.IgnoreCommandLineOptions = True

# Processes to include in the plots
PROCESSES = [
    'DYJets', 'GGToZZ_BSI', 'ST', 'TT', 'TTV', 'WJets',
    'WW', 'WZ', 'ZZ', 'VVV']


def combine_templates(templates_file, channel, syst=None):
    """Construct combined template for all considered processes.

    Take into account requested systematic variation.  It it does not
    affect a process, take the nominal template for it.

    Arguments:
        templates_file:  ROOT file with templates.
        channel:         Label of requested channel.
        syst:            Fully qualified systematic variation.  If None,
            take nominal templates for all processes.

    Return value:
        Hist1D representing the combined template.
    """

    total = Hist1D()
    for process in PROCESSES:
        hist = None
        if syst:
            hist = templates_file.Get('/'.join([channel, process, syst]))
        if not hist:
            hist = templates_file.Get('/'.join([channel, process, 'nominal']))
        total += hist
    return total


def plot(nominal, up, down, save_path, title='', channel_label='',
         formats=['pdf']):
    """Plot impact of a systematic variation.

    Arguments:
        nominal, up, down:  Nominal and alternative templates
            represented with Hist1D.
        save_path:  Path to save produced figure, without file
            extension.
        title:          Title for the figure.
        channel_label:  Label of the channel to be included in the plot.
        formats:        Formats in which to save the figure.

    Return value:
        None.
    """

    fig = plt.figure(figsize=(6.4, 6.4))
    fig.patch.set_alpha(0.)
    gs = mpl.gridspec.GridSpec(2, 1, hspace=0., height_ratios=[1, 1])
    axes_upper = fig.add_subplot(gs[0, 0])
    axes_lower = fig.add_subplot(gs[1, 0])

    binning = nominal.binning
    centres = (binning[:-1] + binning[1:]) / 2
    widths = binning[1:] - binning[:-1]

    # In the upper panel, plot the nominal distribution and
    # uncertainty bands, in the inverse order
    bands = [
        [
            (nominal.contents - nominal.errors)[1:-1],
            (nominal.contents + nominal.errors)[1:-1],
            '0.85'
        ],
        [
            down.contents[1:-1],
            up.contents[1:-1],
            mpl.colors.to_rgba('C0', 0.3)
        ]
    ]
    for band in bands:
        band[0] = band[0] / widths
        band[1] = band[1] / widths
    for band in bands:
        axes_upper.fill_between(
            binning, list(band[0]) + [band[0][-1]],
            list(band[1]) + [band[1][-1]],
            step='post', color=band[2]
        )

    axes_upper.hist(
        binning[:-1], bins=binning,
        weights=nominal.contents[1:-1] / widths,
        histtype='step', color='black'
    )

    # Plot relative deviation in the lower panel
    for hist, direction, marker in [(up, 'up', '^'), (down, 'down', 'v')]:
        deviation = hist.contents[1:-1] / \
            nominal.contents[1:-1] - 1
        axes_lower.plot(
            centres, deviation * 100,
            marker=marker, mfc='C0',
            c=mpl.colors.to_rgba('C0', 0.4),
            label=direction
        )

    for axes in [axes_upper, axes_lower]:
        axes.set_xlim(binning[0], binning[-1])
    axes_upper.set_yscale('log')

    # Restrict the range for the deviation in order to avoid
    # exponentials in the y axis of the lower pad
    ylim = list(axes_lower.get_ylim())
    if ylim[0] > -0.15:
        ylim[0] = -0.15
    if ylim[1] < 0.15:
        ylim[1] = 0.15
    axes_lower.set_ylim(ylim)

    axes_upper.set_title(title)
    axes_upper.set_ylabel(
        r'$\langle\mathrm{Events}\//\/\mathrm{GeV}\rangle$'
    )
    axes_upper.text(
        0., 1.002, channel_label, ha='left', va='bottom',
        transform=axes_upper.transAxes
    )

    # Remove tick labels in the upper axes
    axes_upper.set_xticklabels([])

    axes_lower.axhline(0., ls='dashed', lw=0.8, c='black')
    axes_lower.legend()

    axes_lower.set_xlabel('$m_\\mathrm{T}$ [GeV]')
    axes_lower.set_ylabel('Deviation [%]')

    for fmt in formats:
        fig.savefig(f'{save_path}.{fmt}')
    plt.close(fig)


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(__doc__)
    arg_parser.add_argument('templates', help='ROOT file with templates.')
    arg_parser.add_argument('-o', '--output', default='fig',
                            help='Directory for plots.')
    arg_parser.add_argument('-f', '--formats', nargs='+', default=['pdf'],
                            help='Formats for plots.')
    args = arg_parser.parse_args()

    try:
        os.makedirs(args.output)
    except FileExistsError:
        pass
    plt.style.use(mpl_style)

    templates_file = ROOT.TFile(args.templates)
    for channel, channel_label in [
        ('eq0jets', '$0j$'), ('geq1jets', '$\\geqslant 1j$'), ('vbf', 'VBF')
    ]:
        # Find systematic uncertainties affecting this channel
        systs = set()
        for process in PROCESSES:
            d = templates_file.Get('/'.join([channel, process]))
            for key in d.GetListOfKeys():
                template_name = key.GetName()
                if template_name.endswith('Up'):
                    systs.add(template_name[:-2])

        total_nominal = combine_templates(templates_file, channel)
        for syst in systs:
            total_up = combine_templates(templates_file, channel, syst + 'Up')
            total_down = combine_templates(templates_file, channel,
                                           syst + 'Down')
            plot(
                total_nominal, total_up, total_down,
                os.path.join(args.output, f'{channel}_{syst}'),
                formats=args.formats, title=syst, channel_label=channel_label)

    templates_file.Close()
