use crate::types::PositionalData;

/// Calculates Reynolds' distance (1983): D_R = −ln(1 − GST).
///
/// Linear with drift time for recently diverged populations.
/// Not affected by `--normalize`.
pub fn calculate_reynolds_distance_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &[usize],
) -> f64 {
    let gst = super::gst::calculate_gst_for_pair(data1, data2, all_positions);

    if (1.0 - gst) > 0.0 {
        -(1.0 - gst).ln()
    } else {
        f64::INFINITY
    }
}
