# Changelog

All notable changes to st-speckmol will be documented in this file.

## [0.0.7] - 2024-12-17

### Fixed
- **Major Fix**: Resolved Streamlit Cloud rendering issue ([#22](https://github.com/avrabyt/st-speckmol/issues/22))
  - Molecules now display correctly on Streamlit Cloud
  - Switched from `ipywidgets.embed` to native `ipyspeck.stspeck` Streamlit component
  - Maintains 100% backward compatibility with existing API
  - All parameters and functionality preserved

### Changed
- Updated implementation to use `ipyspeck`'s native Streamlit component when available
- Falls back to ipywidgets approach for local development if needed
- No breaking changes - all existing code continues to work

### Technical Details
- The fix leverages `ipyspeck` v0.6.1's built-in Streamlit component (`stspeck`)
- Uses proper `components.declare_component()` architecture
- Eliminates JavaScript module loading conflicts on Streamlit Cloud

## [0.0.6] - Previous Release
- Added `speck_plot` function
- Deprecated `spec_plot` function
- Parameter handling improvements

## Earlier Versions
See git history for details.

