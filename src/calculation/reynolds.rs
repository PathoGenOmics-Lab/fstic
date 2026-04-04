use crate::types::{GenomicPos, PositionalData};

/// Calculates Reynolds' distance (1983): D_R = −ln(1 − GST).
///
/// Not affected by `--normalize`.
pub fn calculate_reynolds_distance_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &[GenomicPos],
) -> f64 {
    let gst = super::gst::calculate_gst_for_pair(data1, data2, all_positions);
    if (1.0 - gst) > 0.0 { -(1.0 - gst).ln() } else { f64::INFINITY }
}
