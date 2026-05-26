// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Provenance metadata for tabulated optics inputs and generated products.
//!
//! The types in this module are standalone and avoid any dependency on astronomy-
//! specific crates or data pipelines.

use alloc::string::String;

/// Describes where tabulated or derived optical data came from.
///
/// # Examples
///
/// ```rust
/// use optica::data::DataSource;
///
/// let source = DataSource::Computed {
///     name: String::from("rayleigh-fit"),
/// };
///
/// assert!(matches!(source, DataSource::Computed { .. }));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum DataSource {
    /// A literature citation identified by a bibliography key and optional DOI.
    LiteratureCitation {
        /// Local bibliography key or citation identifier.
        bibkey: String,
        /// DOI string when available.
        doi: Option<String>,
    },
    /// A file bundled with the crate or application.
    BundledFile {
        /// Relative or virtual path to the bundled file.
        path: String,
    },
    /// An externally hosted resource.
    External {
        /// URL for the external resource.
        url: String,
    },
    /// Data produced algorithmically at build time or runtime.
    Computed {
        /// Human-readable name of the generating procedure.
        name: String,
    },
}

/// Provenance metadata attached to a table, spectrum, or optical product.
///
/// # Examples
///
/// ```rust
/// use optica::data::Provenance;
///
/// let provenance = Provenance::cited("bodhaine1999")
///     .with_version("v1")
///     .with_notes("Derived from published coefficients.");
///
/// assert_eq!(provenance.version.as_deref(), Some("v1"));
/// assert!(provenance.notes.as_deref().unwrap().contains("published"));
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Provenance {
    /// Origin of the underlying data, if known.
    pub source: Option<DataSource>,
    /// Version string for the source data or generating process.
    pub version: Option<String>,
    /// Retrieval timestamp or release date in caller-defined text form.
    pub retrieved_at: Option<String>,
    /// Optional SHA-256 or similar 32-byte checksum.
    pub checksum: Option<[u8; 32]>,
    /// Free-form notes about derivation, validation, or caveats.
    pub notes: Option<String>,
}

impl Provenance {
    /// Creates empty provenance metadata.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::data::Provenance;
    ///
    /// let provenance = Provenance::new();
    /// assert!(provenance.source.is_none());
    /// ```
    #[must_use]
    pub const fn new() -> Self {
        Self {
            source: None,
            version: None,
            retrieved_at: None,
            checksum: None,
            notes: None,
        }
    }

    /// Creates provenance for a bundled file.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::data::{DataSource, Provenance};
    ///
    /// let provenance = Provenance::bundled_file("tables/rayleigh.csv");
    /// assert!(matches!(
    ///     provenance.source,
    ///     Some(DataSource::BundledFile { .. })
    /// ));
    /// ```
    #[must_use]
    pub fn bundled_file(path: impl Into<String>) -> Self {
        Self {
            source: Some(DataSource::BundledFile { path: path.into() }),
            ..Self::new()
        }
    }

    /// Creates provenance for a literature citation.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::data::{DataSource, Provenance};
    ///
    /// let provenance = Provenance::cited("bodhaine1999");
    /// assert!(matches!(
    ///     provenance.source,
    ///     Some(DataSource::LiteratureCitation { .. })
    /// ));
    /// ```
    #[must_use]
    pub fn cited(bibkey: impl Into<String>) -> Self {
        Self {
            source: Some(DataSource::LiteratureCitation {
                bibkey: bibkey.into(),
                doi: None,
            }),
            ..Self::new()
        }
    }

    /// Creates provenance for computed data.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::data::{DataSource, Provenance};
    ///
    /// let provenance = Provenance::computed("analytic fit");
    /// assert!(matches!(
    ///     provenance.source,
    ///     Some(DataSource::Computed { .. })
    /// ));
    /// ```
    #[must_use]
    pub fn computed(name: impl Into<String>) -> Self {
        Self {
            source: Some(DataSource::Computed { name: name.into() }),
            ..Self::new()
        }
    }

    /// Attaches a version string.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::data::Provenance;
    ///
    /// let provenance = Provenance::new().with_version("2026-01");
    /// assert_eq!(provenance.version.as_deref(), Some("2026-01"));
    /// ```
    #[must_use]
    pub fn with_version(mut self, version: impl Into<String>) -> Self {
        self.version = Some(version.into());
        self
    }

    /// Attaches free-form notes.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::data::Provenance;
    ///
    /// let provenance = Provenance::new().with_notes("Interpolated onto a uniform grid.");
    /// assert!(provenance.notes.as_deref().unwrap().contains("uniform"));
    /// ```
    #[must_use]
    pub fn with_notes(mut self, notes: impl Into<String>) -> Self {
        self.notes = Some(notes.into());
        self
    }
}

/// Identifies how a table's numeric payload is stored or generated.
///
/// # Examples
///
/// ```rust
/// use optica::data::TableSource;
///
/// static VALUES: [f64; 2] = [1.0, 2.0];
/// let source = TableSource::EmbeddedSlices {
///     data: &VALUES,
///     len: VALUES.len(),
/// };
///
/// assert!(matches!(source, TableSource::EmbeddedSlices { len: 2, .. }));
/// ```
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub enum TableSource {
    /// Opaque embedded bytes, suitable for a custom parser.
    EmbeddedBytes(&'static [u8]),
    /// Pre-decoded embedded floating-point samples.
    EmbeddedSlices {
        /// Numeric data buffer.
        data: &'static [f64],
        /// Logical length of the encoded table.
        len: usize,
    },
    /// Path to an external table on disk.
    ExternalPath(String),
    /// A generated table described by text.
    Generated {
        /// Human-readable description of the generator.
        description: String,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn provenance_builders_compose() {
        let provenance = Provenance::bundled_file("data/fit.csv")
            .with_version("v1")
            .with_notes("validated");

        assert_eq!(provenance.version.as_deref(), Some("v1"));
        assert_eq!(provenance.notes.as_deref(), Some("validated"));
    }
}
