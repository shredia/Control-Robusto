"""Configuración de alto nivel para los experimentos de control.

Este módulo es la única fuente que decide qué estimador y compensadores se usan.
Los parámetros físicos y ganancias numéricas permanecen en ``parameters.py``.
"""

from dataclasses import dataclass, replace
from typing import ClassVar


@dataclass(frozen=True)
class ProjectConfig:
    """Selección de arquitectura para una simulación.

    ``observer`` admite:
      - ``real``: sensores ideales, útil como referencia.
      - ``ekf``: EKF electromecánico sin corrección angular HFI.
      - ``ekf_hfi``: EKF electromecánico con corrección HFI.
      - ``papers``: ángulo HFI + Kalman cinemático [theta, omega, alpha].
    """

    observer: str = "real"
    enable_dob: bool = False
    enable_lms: bool = False
    enable_lob: bool = False
    enable_bpf: bool = True
    sensorless_control: bool = False
    inverter_mode: str = "average"
    output_name: str = "simulation"

    OBSERVERS: ClassVar[tuple[str, ...]] = ("real", "ekf", "ekf_hfi", "papers")

    def validate(self) -> "ProjectConfig":
        if self.observer not in self.OBSERVERS:
            raise ValueError(
                f"observer={self.observer!r} no es válido; use {', '.join(self.OBSERVERS)}"
            )
        if self.inverter_mode not in ("average", "switching"):
            raise ValueError("inverter_mode debe ser 'average' o 'switching'")
        if self.enable_lms and not self.enable_dob:
            raise ValueError("LMS necesita DOB activado porque aprende desde su perturbación estimada")
        if self.sensorless_control and self.observer == "real":
            raise ValueError("sensorless_control no puede usar observer='real'")
        return self

    @property
    def enable_hfi(self) -> bool:
        return self.observer in ("ekf_hfi", "papers")

    @property
    def feedback_source(self) -> str:
        if not self.sensorless_control:
            return "real"
        return "papers" if self.observer == "papers" else "ekf"

    def with_changes(self, **changes) -> "ProjectConfig":
        return replace(self, **changes).validate()


PRESETS = {
    "sensed": ProjectConfig(output_name="sensed"),
    "ekf": ProjectConfig(observer="ekf", sensorless_control=True, output_name="ekf"),
    "ekf-dob": ProjectConfig(
        observer="ekf", enable_dob=True, sensorless_control=True, output_name="ekf_dob"
    ),
    "ekf-hfi": ProjectConfig(
        observer="ekf_hfi", sensorless_control=True, output_name="ekf_hfi"
    ),
    "papers": ProjectConfig(
        observer="papers", sensorless_control=True, output_name="papers"
    ),
    "papers-dob-lms": ProjectConfig(
        observer="papers", enable_dob=True, enable_lms=True,
        sensorless_control=True, output_name="papers_dob_lms"
    ),
}


def get_preset(name: str) -> ProjectConfig:
    try:
        return PRESETS[name].validate()
    except KeyError as exc:
        raise ValueError(f"Perfil desconocido {name!r}; use {', '.join(PRESETS)}") from exc
