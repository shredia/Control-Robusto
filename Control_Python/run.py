"""Entrada única para ejecutar y comparar arquitecturas del proyecto."""

import argparse
from dataclasses import replace
from pathlib import Path

import numpy as np

from project_config import PRESETS, ProjectConfig, get_preset
from simulation import PMSMSimulation


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "outputs"


def position_reference(t):
    return 10.0 * t


def speed_reference(t):
    return 10.0


def load_torque(t, theta_m=0.0):
    step = 0.02 if t >= 0.5 else 0.0
    gravity = 0.01 * np.sin(theta_m)
    return step + gravity


def build_parser():
    parser = argparse.ArgumentParser(
        description="Simulación configurable del motor y sus observadores."
    )
    parser.add_argument(
        "--profile", choices=tuple(PRESETS), default="sensed",
        help="Arquitectura base (por defecto: sensed).",
    )
    parser.add_argument(
        "--observer", choices=ProjectConfig.OBSERVERS,
        help="Sobrescribe el observador del perfil.",
    )
    parser.add_argument("--dob", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--lms", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--lob", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--bpf", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument(
        "--sensorless", action=argparse.BooleanOptionalAction, default=None,
        help="Cierra los lazos con el observador seleccionado.",
    )
    parser.add_argument("--duration", type=float, help="Duración de simulación en segundos.")
    parser.add_argument("--plots", action="store_true", help="Guarda los gráficos PNG.")
    return parser


def resolve_config(args):
    config = get_preset(args.profile)
    changes = {}
    mappings = {
        "observer": args.observer,
        "enable_dob": args.dob,
        "enable_lms": args.lms,
        "enable_lob": args.lob,
        "enable_bpf": args.bpf,
        "sensorless_control": args.sensorless,
    }
    for field, value in mappings.items():
        if value is not None:
            changes[field] = value
    return replace(config, **changes).validate()


def run(config, duration=None):
    print("\nConfiguración activa")
    print(f"  Observador : {config.observer}")
    print(f"  Sensorless : {config.sensorless_control}")
    print(f"  HFI / BPF  : {config.enable_hfi} / {config.enable_bpf}")
    print(f"  DOB / LMS  : {config.enable_dob} / {config.enable_lms}")
    print(f"  LOB        : {config.enable_lob}\n")

    sim = PMSMSimulation(t_end=duration, config=config)
    return sim.run(position_reference, load_torque, speed_reference)


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        config = resolve_config(args)
    except ValueError as exc:
        raise SystemExit(f"Configuración inválida: {exc}") from exc

    history = run(config, duration=args.duration)
    if args.plots:
        # Se reutilizan temporalmente los gráficos del script histórico.
        from main import plot_bpf_zoom, plot_results
        plot_bpf_zoom(history, title_prefix=config.output_name)
        plot_results(history, title_prefix=config.output_name)

    print("\nResumen")
    print(f"  RMSE velocidad ({config.observer}): {history['rmse_w_observer']:.4f} rad/s")
    print(f"  RMSE ángulo seleccionado          : {history['rmse_theta_observer']:.4f} rad_e")
    return history


if __name__ == "__main__":
    main()
