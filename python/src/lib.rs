use pyo3::prelude::*;
use mynum::mpz::Mpz;
use mynum::mpf::Mpf;

#[pyclass(name = "Mpz")]
#[derive(Clone)]
struct PyMpz {
    inner: Mpz,
}

#[pymethods]
impl PyMpz {
    #[new]
    fn new(s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: Mpz::from_str(s, 10)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?,
        })
    }

    fn __add__(&self, other: &PyMpz) -> PyMpz {
        PyMpz { inner: self.inner.add(&other.inner) }
    }

    fn __sub__(&self, other: &PyMpz) -> PyMpz {
        PyMpz { inner: self.inner.sub(&other.inner) }
    }

    fn __mul__(&self, other: &PyMpz) -> PyMpz {
        PyMpz { inner: self.inner.mul(&other.inner) }
    }

    fn __str__(&self) -> String {
        self.inner.to_string(10)
    }

    fn __repr__(&self) -> String {
        format!("Mpz('{}')", self.inner.to_string(10))
    }

    fn __eq__(&self, other: &PyMpz) -> bool {
        self.inner.cmp(&other.inner) == core::cmp::Ordering::Equal
    }
}

#[pyclass(name = "Mpf")]
#[derive(Clone)]
struct PyMpf {
    inner: Mpf,
}

#[pymethods]
impl PyMpf {
    #[new]
    fn new(value: f64, precision: Option<usize>) -> Self {
        Self { inner: Mpf::from_f64(value, precision.unwrap_or(64)) }
    }

    #[staticmethod]
    fn from_string(s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: Mpf::from_str(s, 10)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?,
        })
    }

    fn __add__(&self, other: &PyMpf) -> PyMpf {
        PyMpf { inner: self.inner.add(&other.inner) }
    }

    fn __mul__(&self, other: &PyMpf) -> PyMpf {
        PyMpf { inner: self.inner.mul(&other.inner) }
    }

    fn __str__(&self) -> String {
        self.inner.to_string(10)
    }

    fn __repr__(&self) -> String {
        format!("Mpf('{}')", self.inner.to_string(10))
    }
}

#[pymodule]
fn mynum_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMpz>()?;
    m.add_class::<PyMpf>()?;
    Ok(())
}
