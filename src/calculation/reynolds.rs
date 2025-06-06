use crate::types::PositionalData;
use std::collections::HashSet;

/// Calculates Reynolds' distance (1983), which is a transformation of GST.
/// D = -ln(1 - GST). This distance is theoretically linear with time.
/// It is a global metric and not affected by the --normalize flag.
pub fn calculate_reynolds_distance_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &HashSet<usize>,
) -> f64 {
    // This distance is a direct function of the global GST.
    // To avoid code duplication, we call the GST calculation function directly.
    let gst = super::gst::calculate_gst_for_pair(data1, data2, all_positions);

    if (1.0 - gst) > 0.0 {
        -(1.0 - gst).ln()
    } else {
        f64::INFINITY // Distance is infinite if populations are completely differentiated (GST=1)
    }
}
