"""命令行读取一个具名设计配置并导出计算结果。"""
import argparse
from pathlib import Path

from .api import analyze_design
from .io import read_config, export_result, write_json


def main(argv=None):
    parser = argparse.ArgumentParser(description="直齿／弧齿锥齿轮几何、刀具、强度与三维齿面")
    parser.add_argument("config", help="JSON设计配置；cad-import模式为STEP文件")
    parser.add_argument("--output", default="output/design", help="结果目录")
    parser.add_argument("--plots", action="store_true", help="另外导出PNG/SVG工程图")
    parser.add_argument("--view3d", action="store_true", help="另外导出离线三维装配HTML和场景JSON，无需绘图库")
    parser.add_argument("--surface3d", action="store_true", help="导出理想直齿解析齿面查看页及OBJ/JSON，不含根过渡或实体轮体")
    parser.add_argument("--cad",action="store_true",help="展成模式用可选OCP内核导出五面壳STEP/BREP及三角网格")
    parser.add_argument("--render3d",action="store_true",help="展成CAD模式用可选VTK内核导出三维PNG/VTP")
    parser.add_argument("--publication",action="store_true",help="图表另外导出论文用PDF，保持单/双栏物理尺寸")
    parser.add_argument("--mode", choices=("design","manufacturing","generated","cad-import"),default="design",help="输入配置所属功能，默认design")
    arguments = parser.parse_args(argv)
    if arguments.mode!="design" and (arguments.view3d or arguments.surface3d):
        parser.error("--view3d/--surface3d只用于零偏置design配置")
    if arguments.plots and arguments.mode in ("manufacturing","cad-import"):
        parser.error("图表选项用于design或generated；其他模式导出具名数据")
    if arguments.cad and arguments.mode not in ('generated','cad-import'):parser.error('--cad用于generated或cad-import模式')
    if arguments.render3d and not(arguments.mode=='cad-import' or (arguments.cad and arguments.mode=='generated')):parser.error('--render3d用于cad-import或同时指定--cad的generated模式')
    if arguments.publication and not arguments.plots:parser.error('--publication须与--plots一起使用')
    try:
        config={} if arguments.mode=='cad-import' else read_config(arguments.config)
        if arguments.mode=="design":
            result = analyze_design(config)
            export_result(result, arguments.output)
            if arguments.view3d or arguments.surface3d:
                from .assembly import export_assembly_view
                export_assembly_view(result, arguments.output,with_tooth_surfaces=arguments.surface3d)
            if arguments.plots:
                from .plotting import plot_design
                plot_design(result, arguments.output,publication=arguments.publication)
        elif arguments.mode=="manufacturing":
            from .manufacturing import calculate_manufacturing,export_manufacturing_result
            result=calculate_manufacturing(config.get("manufacturing",config))
            export_manufacturing_result(result,arguments.output)
        elif arguments.mode=='generated':
            from .generated_surface import build_generated_surface_grid,build_generated_accuracy_grid
            from .generated_io import export_generated_grid
            request=config.get('generated_geometry',config)
            result=build_generated_surface_grid(request,**config.get('sampling',{}))
            export_generated_grid(result,arguments.output)
            if arguments.plots:
                from .plotting import plot_generated_geometry
                plot_generated_geometry(result,arguments.output,publication=arguments.publication)
            if arguments.cad:
                from .cad import build_generated_shell,export_step,export_brep,triangulate,measure_surface_accuracy
                shell=build_generated_shell(result)
                export_step(shell,Path(arguments.output)/'generated_shell.step')
                export_brep(shell,Path(arguments.output)/'generated_shell.brep')
                mesh=triangulate(shell)
                write_json(mesh,Path(arguments.output)/'generated_cad_mesh.json')
                write_json(dict(shell=shell.diagnostics,accuracy=measure_surface_accuracy(shell,build_generated_accuracy_grid(request))),
                           Path(arguments.output)/'generated_cad_diagnostics.json')
                if arguments.render3d:
                    from .cad_preview import render_cad_mesh,export_vtp
                    render_cad_mesh(mesh,Path(arguments.output)/'generated_cad.png')
                    export_vtp(mesh,Path(arguments.output)/'generated_cad.vtp')
        elif arguments.mode=='cad-import':
            from .cad_import import load_step_model,triangulate_imported_model,export_model_step,export_model_brep
            model=load_step_model(arguments.config)
            result=model.diagnostics
            write_json(result,Path(arguments.output)/'imported_cad_diagnostics.json')
            mesh=triangulate_imported_model(model)
            write_json(mesh,Path(arguments.output)/'imported_cad_mesh.json')
            if arguments.cad:
                export_model_step(model,Path(arguments.output)/'imported_model.step')
                export_model_brep(model,Path(arguments.output)/'imported_model.brep')
            if arguments.render3d:
                from .cad_preview import render_cad_mesh,export_vtp
                render_cad_mesh(mesh,Path(arguments.output)/'imported_cad.png')
                export_vtp(mesh,Path(arguments.output)/'imported_cad.vtp')
    except ModuleNotFoundError as error:
        if error.name in ("numpy", "matplotlib"):
            parser.exit(2, '绘图需要可选依赖，请先运行 python -m pip install ".[plots]"\n')
        if error.name in ('OCP','vtk','vtkmodules') or str(error.name).startswith(('OCP.','vtkmodules.')):
            parser.exit(2, 'CAD/三维渲染需要可选依赖，请先运行 python -m pip install ".[cad,render3d]"\n')
        raise
    except ImportError as error:
        parser.exit(2,f"可选几何/渲染依赖不可用：{error}\n请按需安装 .[cad,render3d]\n")
    except (OSError, ValueError, KeyError, TypeError, ArithmeticError, NotImplementedError) as error:
        parser.exit(2, f"输入输出或模型范围错误：{error}\n")
    if arguments.mode=="design":
        rating = result["rating"]
        print(f"接触应力：{rating['contact_stress_pa']/1e6:.3f} MPa")
        print(f"两轮弯曲应力：{rating['pinion_bending_stress_pa']/1e6:.3f} / {rating['wheel_bending_stress_pa']/1e6:.3f} MPa")
        print(f"满足本次安全系数目标：{'是' if rating['meets_targets'] else '否'}")
    else:
        print(f"功能计算完成：{arguments.mode}")
    print(f"结果已保存至：{arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
