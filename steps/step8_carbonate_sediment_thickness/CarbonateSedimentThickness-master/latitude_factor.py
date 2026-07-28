"""
Latitude- and time-dependent weighting for pelagic carbonate accumulation
rate.

This module provides ``latitude_factor()``, a multiplier in [0, 1] that is
applied to the depth-dependent carbonate sedimentation rate computed in
``carbonate_sediment_thickness.calc_carbonate_decompacted_sediment_thickness``.
It corrects the spurious accumulation of thick pelagic carbonate on
young, shallow, high-latitude crust (most visibly the Arctic Ocean) that
arises because the underlying algorithm has no temperature, productivity
or carbonate-saturation-state control - only bathymetry relative to a
globally uniform CCD.

The correction is *time-dependent*. The strong latitudinal cut-off of
pelagic carbonate is an icehouse phenomenon: it emerges only as Antarctic
glaciation begins (Eocene-Oligocene Transition, ~34 Ma) and reaches its
modern form during the Oligocene CO2 drawdown (~34 -> ~23 Ma, Hoenisch
et al. 2023). Before ~34 Ma the high southern latitudes were producing
pelagic carbonate (Maud Rise ODP 689/690, Mentelle Basin IODP U1513-U1516,
Campbell Plateau DSDP 277), so no latitudinal restriction is applied to
times >= 34 Ma. Between 34 Ma and 23 Ma the strength of the restriction
ramps smoothly in via a half-cosine in time. At <= 23 Ma the full
modern latitudinal restriction is applied.

References (modern latitudinal pattern - applied at t <= 23 Ma):
    Archer, D. (1996a). An atlas of the distribution of calcium carbonate
        in sediments of the deep sea. Global Biogeochemical Cycles, 10(1),
        159-174. doi:10.1029/95GB03016
    Archer, D. (1996b). A data-driven model of the calcite lysocline.
        Global Biogeochemical Cycles, 10(3), 511-526. doi:10.1029/96GB01521
    Honjo, S., Manganini, S. J., Krishfield, R. A., & Francois, R. (2008).
        Particulate organic carbon fluxes to the ocean interior and factors
        controlling the biological pump: A synthesis of global sediment
        trap programs since 1983. Progress in Oceanography, 76(3), 217-285.
        doi:10.1016/j.pocean.2007.11.003
    Jutterstrom, S., & Anderson, L. G. (2016). Disparate acidification and
        calcium carbonate desaturation of deep and shallow waters of the
        Arctic Ocean. Nature Communications, 7, 12821.
        doi:10.1038/ncomms12821
    Schiebel, R. (2002). Planktic foraminiferal sedimentation and the
        marine calcite budget. Global Biogeochemical Cycles, 16(4), 1065.
        doi:10.1029/2001GB001459
    Stein, R. (2008). Arctic Ocean Sediments: Processes, Proxies, and
        Paleoenvironment. Developments in Marine Geology, Vol. 2, Elsevier.
    Sulpis, O., Boudreau, B. P., Mucci, A., Jenkins, C., Trossman, D. S.,
        Arbic, B. K., & Key, R. M. (2018). Current CaCO3 dissolution at the
        seafloor caused by anthropogenic CO2. PNAS, 115(46), 11700-11705.
        doi:10.1073/pnas.1804250115

References (warm-climate high-latitude pelagic carbonate - no restriction
at t >= 34 Ma):
    Thomas, E. (1990). Late Cretaceous through Neogene deep-sea benthic
        foraminifers (Maud Rise, Weddell Sea, Antarctica). Proceedings of
        the Ocean Drilling Program, Scientific Results, 113, 571-594.
        (ODP Sites 689 and 690: Maastrichtian-Cenozoic pelagic carbonate at
        ~65 deg S palaeolatitude.)
    Huber, B. T. (1991). Maestrichtian planktonic foraminifer biostratigraphy
        and the Cretaceous/Tertiary boundary at ODP Hole 738C (Kerguelen
        Plateau, southern Indian Ocean). Proceedings of the Ocean Drilling
        Program, Scientific Results, 119, 451-465.
    Huber, B. T., Hobbs, R. W., Bogus, K. A., & the Expedition 369
        Scientists (2019). Australia Cretaceous Climate and Tectonics.
        Proceedings of the IODP, Volume 369. doi:10.14379/iodp.proc.369.2019
        (IODP Sites U1513-U1516, Mentelle Basin: pelagic carbonate at
        54-60 deg S palaeolatitude through OAE 2 and the PETM.)
    Gallagher, S. J., Wagstaff, B. E., Baird, J. G., Wallace, M. W., &
        Li, C. L. (2020). Eocene to Oligocene high palaeolatitude neritic
        record of Oi-1 glaciation in the Otway Basin, southeast Australia.
        Global and Planetary Change, 191, 103217.
        doi:10.1016/j.gloplacha.2020.103217

References (timing of Antarctic glaciation onset and CO2 drawdown):
    Coxall, H. K., Wilson, P. A., Palike, H., Lear, C. H., & Backman, J.
        (2005). Rapid stepwise onset of Antarctic glaciation and deeper
        calcite compensation in the Pacific Ocean. Nature, 433, 53-57.
        doi:10.1038/nature03135
    Miller, K. G., Wright, J. D., Katz, M. E., Browning, J. V.,
        Cramer, B. S., Wade, B. S., & Mizintseva, S. F. (2008). A view of
        Antarctic ice-sheet evolution from sea-level and deep-sea isotope
        changes during the Late Cretaceous-Cenozoic. In Antarctica: A
        Keystone in a Changing World (pp. 55-70). NAS Press.
    Taylor, V. E., Wilson, P. A., Bohaty, S. M., Foster, G. L., Lear, C. H.,
        & Cooper, M. J. (2023). Transient shoaling, over-deepening and
        settling of the calcite compensation depth at the Eocene-Oligocene
        transition. Paleoceanography and Paleoclimatology, 38, e2022PA004493.
        doi:10.1029/2022PA004493
    CenCO2PIP Consortium, Hoenisch, B., et al. (2023). Toward a Cenozoic
        history of atmospheric CO2. Science, 382, eadi5177.
        doi:10.1126/science.adi5177
        (Constrains the ~34 -> ~23 Ma pCO2 decline that motivates the
        temporal ramp used here.)
"""

import math


# Latitude (degrees, absolute value) at and equatorward of which no
# high-latitude reduction is applied. Equivalent to the poleward edge of the
# productive subpolar belt in the modern ocean.
HIGH_LATITUDE_FULL_PRESERVATION_DEG = 60.0

# Latitude (degrees, absolute value) at and poleward of which pelagic
# carbonate accumulation is set to zero. Modern deep Arctic basins are
# essentially carbonate-free (Stein 2008) and deep waters are at/below
# calcite saturation (Jutterstrom & Anderson 2016).
HIGH_LATITUDE_ZERO_ACCUMULATION_DEG = 80.0

# Onset of Antarctic glaciation (Eocene-Oligocene Transition, EOT-1 / Oi-1):
# at and before this age no latitudinal restriction is applied. Pelagic
# carbonate was being deposited at >55 deg S palaeolatitude (Maud Rise,
# Mentelle Basin, Campbell Plateau) through the Mesozoic and Paleogene
# greenhouse, so a modern-calibrated restriction would over-strip the
# pre-icehouse record.
GLACIATION_ONSET_MA = 34.0

# End of the Oligocene CO2 drawdown (Oligocene-Miocene boundary):
# at and after this age the full modern latitudinal restriction is applied.
# Hoenisch et al. (2023) constrain pCO2 to have fallen from ~900 ppmv at
# the EOT to ~300-400 ppmv by ~23 Ma, by which point the climate state is
# recognisably modern and the latitudinal cut-off is fully expressed.
FULL_CORRECTION_MA = 23.0


def _spatial_factor(lat_deg, lat_full, lat_zero):
    """Half-cosine taper in absolute latitude. Returns 1 equatorward of
    ``lat_full``, 0 poleward of ``lat_zero``, smooth in between.
    """
    abs_lat = abs(lat_deg)
    if abs_lat <= lat_full:
        return 1.0
    if abs_lat >= lat_zero:
        return 0.0
    x = (abs_lat - lat_full) / (lat_zero - lat_full)
    return 0.5 * (1.0 + math.cos(math.pi * x))


def _temporal_blend(time_ma, glaciation_onset_ma, full_correction_ma):
    """Half-cosine ramp controlling how strongly the spatial restriction
    is applied as a function of paleo-time.

    Returns 0 at and before ``glaciation_onset_ma`` (greenhouse - no
    restriction), 1 at and after ``full_correction_ma`` (icehouse - full
    modern restriction), smooth in between.
    """
    if time_ma >= glaciation_onset_ma:
        return 0.0
    if time_ma <= full_correction_ma:
        return 1.0
    x = ((glaciation_onset_ma - time_ma)
         / (glaciation_onset_ma - full_correction_ma))
    return 0.5 * (1.0 - math.cos(math.pi * x))


def latitude_factor(
        lat_deg,
        time_ma=0.0,
        lat_full=HIGH_LATITUDE_FULL_PRESERVATION_DEG,
        lat_zero=HIGH_LATITUDE_ZERO_ACCUMULATION_DEG,
        glaciation_onset_ma=GLACIATION_ONSET_MA,
        full_correction_ma=FULL_CORRECTION_MA):
    """
    Latitude- and time-dependent weighting for pelagic carbonate
    accumulation rate.

    Returns a factor in [0, 1] that should multiply the depth-dependent
    carbonate sedimentation rate. Formed by combining:
      (a) a spatial half-cosine taper in absolute latitude that is 1
          equatorward of ``lat_full``, 0 poleward of ``lat_zero``, and
      (b) a temporal half-cosine ramp that switches the spatial taper off
          at and before ``glaciation_onset_ma`` (greenhouse), and fully on
          at and after ``full_correction_ma`` (icehouse).

    The combined factor is
        f(|phi|, t) = 1 - b(t) * (1 - f_lat(|phi|))
    which equals 1 everywhere when b(t) = 0 (no effect on greenhouse
    record), equals f_lat(|phi|) when b(t) = 1 (modern restriction), and
    is a smooth blend between the two during the EOT -> Oligocene-Miocene
    transition.

    Parameters
    ----------
    lat_deg : float
        Paleo-latitude (in degrees) of the reconstructed point at
        ``time_ma``. Must be the paleo-latitude rather than the present-day
        latitude, so that crust which originated at low latitude still
        receives credit for carbonate it accumulated before drifting
        poleward.
    time_ma : float, optional
        Paleo-time (Ma) of the reconstructed point. Default 0 (present
        day; full modern restriction). The carbonate workflow calls this
        with ``reconstruction_time`` so the correction is automatically
        switched off in greenhouse intervals.
    lat_full : float, optional
        Absolute latitude (degrees) at and equatorward of which the
        spatial factor is 1. Default 60.
    lat_zero : float, optional
        Absolute latitude (degrees) at and poleward of which the spatial
        factor is 0. Default 80.
    glaciation_onset_ma : float, optional
        Age (Ma) at and before which no latitudinal restriction is
        applied. Default 34.0 (EOT-1 / Oi-1 Antarctic glaciation onset).
    full_correction_ma : float, optional
        Age (Ma) at and after which the full modern restriction is
        applied. Default 23.0 (Oligocene-Miocene boundary, end of the
        Hoenisch et al. 2023 CO2 drawdown).

    Returns
    -------
    float
        Combined weighting factor in [0, 1].

    Rationale
    ---------
    The depth-only formulation in ``carbonate_sediment_thickness`` has no
    temperature, productivity or carbonate-saturation-state control. Without
    a latitude weight it deposits spurious thick carbonate on young,
    shallow high-latitude crust. But the strong latitudinal cut-off of
    pelagic carbonate is an *icehouse* phenomenon, controlled by cold
    surface waters that are weakly saturated with respect to calcite, by
    the regional shoaling of the lysocline/CCD at high latitudes, and
    (in the Arctic specifically) by ice-rafted-debris and clay dilution.
    These conditions only emerge in their modern form during the
    Cenozoic icehouse. In the pre-34 Ma greenhouse, high palaeolatitudes
    were producing pelagic carbonate:
      - Maud Rise (ODP 689/690): Maastrichtian-Eocene calcareous oozes
        at ~65 deg S palaeolatitude (Thomas 1990; Huber 1991).
      - Mentelle Basin (IODP U1513-U1516): continuous pelagic carbonate
        across OAE 2 and the PETM at 54-60 deg S palaeolatitude
        (Huber et al. 2019).
      - Campbell Plateau (DSDP 277) and the Otway Basin (Gallagher et
        al. 2020): carbonate -> non-carbonate turnover at high palaeo-
        latitudes localised exactly across the EOT.

    The temporal ramp is anchored on the EOT-1 / Oi-1 events (Coxall et
    al. 2005; Miller et al. 2008) at ~34 Ma, and on the end of the
    Oligocene pCO2 drawdown at ~23 Ma reconstructed by the CenCO2PIP
    consortium (Hoenisch et al. 2023). A half-cosine ramp (rather than
    linear) is used so the time-derivative of the correction is zero at
    both ends of the window, avoiding artefacts in time-derivative
    diagnostics (e.g. fluxes, plate-influx rates).

    Sensitivity tests
    -----------------
    Spatial defaults (60, 80) are the published recommendation. Temporal
    defaults (34, 23) can be tightened or relaxed, e.g.:
      (34, 33): correction switched on as a step at the Oi-1 event.
      (34, 16): correction ramps over the entire Oligocene-Miocene
                cooling that culminates in the mid-Miocene Climatic
                Optimum / Climate Transition.
    """
    f_lat = _spatial_factor(lat_deg, lat_full, lat_zero)
    b = _temporal_blend(time_ma, glaciation_onset_ma, full_correction_ma)
    return 1.0 - b * (1.0 - f_lat)
