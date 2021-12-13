//
// Created by tomlucas on 06.04.20.
//

#include "PoseRenderer.h"

#include <vtkLine.h>
#include <vtkGlyph3D.h>
#include <vtkSphereSource.h>

#define LOCK std::lock_guard<std::mutex> lockGuard(mt);
namespace adekf::viz
{

    vtkSmartPointer<vtkActor> PoseRenderer::createActor(const char *color, double scale)
    {
        auto colors =
            vtkSmartPointer<vtkNamedColors>::New();
        // Create a sphere
        auto cubeSource =
            vtkSmartPointer<vtkCubeSource>::New();
        cubeSource->SetCenter(0.0, 0.0, 0.0);
        cubeSource->SetBounds(-1 * scale, 1 * scale, -0.5 * scale, 0.5 * scale, -0.25 * scale, 0.25 * scale);
        cubeSource->Update();
        // Create a mapper and actor
        auto mapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(cubeSource->GetOutputPort());
        auto actor =
            vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(colors->GetColor3d(color).GetData());
        renderer->AddActor(actor);
        return actor;
    }

    vtkSmartPointer<vtkActor> PoseRenderer::addPositionBubble(const Eigen::Vector3d &position, const char *color, double radius)
    {
        auto colors =
            vtkSmartPointer<vtkNamedColors>::New();
        // Create a sphere
        auto sphereSource =
            vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(position.x(), position.y(), position.z());
        sphereSource->SetRadius(radius);
        sphereSource->Update();
        // Create a mapper and actor
        auto mapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(sphereSource->GetOutputPort());
        auto actor =
            vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(colors->GetColor3d(color).GetData());
        renderer->AddActor(actor);
        return actor;
    }

    bool PoseRenderer::isDone()
    {
        return window->GetDone();
    }

    void PoseRenderer::setCameraPosition(int x, int y, int z)
    {
        renderer->GetActiveCamera()->SetPosition(x, y, z);
    }

    void PoseRenderer::setFocalPoint(int x, int y, int z)
    {
        renderer->GetActiveCamera()->SetFocalPoint(x, y, z);
    }

    void PoseRenderer::initGui()
    {
        auto colors =
            vtkSmartPointer<vtkNamedColors>::New();
        // Create a renderer, render window, and interactor
        renderer =
            vtkSmartPointer<vtkRenderer>::New();
        auto renderWindow =
            vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        window =
            vtkSmartPointer<vtkRenderWindowInteractor>::New();
        window->SetRenderWindow(renderWindow);
        window->Initialize();
        renderWindow->SetSize(800, 800);
        renderer->ResetCamera();
        renderer->GetActiveCamera()->SetPosition(0, 0, 80);
        renderer->GetActiveCamera()->SetFocalPoint(0, 0, 0);

        // Add the actor to the scene
        renderer->SetBackground(colors->GetColor3d("AliceBlue").GetData());
        renderWindow->MakeCurrent();
        // Render and interact

        vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
            vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview
        window->SetInteractorStyle(style);
        renderWindow->Render();
    }

    void PoseRenderer::updateWindow()
    {
        LOCK for (auto poseActor : posesAndActors)
        {
            auto actor = poseActor.second;
            auto pose = poseActor.first;
            auto position = pose->getPosition();
            actor->SetPosition(position.x(), position.y(), position.z());
            auto q = pose->getOrientation();
            actor->SetOrientation(0, 0, 0);
            Eigen::AngleAxis<decltype(q)::Scalar> aa(q);
            auto axis = aa.axis();
            actor->RotateWXYZ(aa.angle() * 180. / M_PI, axis.x(), axis.y(), axis.z());
        }
        if (focus.get() != NULL)
        {
            auto pos = focus->getPosition();
            setCameraPosition(pos.x() + view_vector.x(), pos.y() + view_vector.y(), pos.z() + view_vector.z());
            setFocalPoint(pos.x(), pos.y(), pos.z());
        }
        window->ProcessEvents();
        window->Render();
    }

    void PoseRenderer::disposeWindow()
    {
        for (auto poseActor : posesAndActors)
        {
            delete poseActor.first;
        }
    }

    void PoseRenderer::displayPath(const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> &path, const char *color)
    {
        vtkSmartPointer<vtkPoints> points =
            vtkSmartPointer<vtkPoints>::New();
        for (Eigen::Vector3d vector : path)
        {
            points->InsertNextPoint(vector.data());
        }
        vtkSmartPointer<vtkCellArray> lines =
            vtkSmartPointer<vtkCellArray>::New();

        for (unsigned int i = 0; i < points->GetNumberOfPoints() - 1; i++)
        {
            vtkSmartPointer<vtkLine> line =
                vtkSmartPointer<vtkLine>::New();
            line->GetPointIds()->SetId(0, i);
            line->GetPointIds()->SetId(1, i + 1);
            lines->InsertNextCell(line);
        }

        // Create a polydata to store everything in
        vtkSmartPointer<vtkPolyData> linesPolyData =
            vtkSmartPointer<vtkPolyData>::New();

        // Add the points to the dataset
        linesPolyData->SetPoints(points);

        // Add the lines to the dataset
        linesPolyData->SetLines(lines);

        // Setup actor and mapper
        vtkSmartPointer<vtkNamedColors> colors =
            vtkSmartPointer<vtkNamedColors>::New();

        vtkSmartPointer<vtkPolyDataMapper> mapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(linesPolyData);

        vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetLineWidth(4);
        actor->GetProperty()->SetColor(colors->GetColor3d(color).GetData());
        renderer->AddActor(actor);
    }

    void PoseRenderer::removeActor(const vtkSmartPointer<vtkActor> &actor)
    {
        LOCK
            renderer->RemoveActor(actor);
    }

    void PoseRenderer::removeActor(const std::pair<PoseReader *, vtkSmartPointer<vtkActor>> &actor)
    {
        LOCK
            posesAndActors.remove(actor);
        renderer->RemoveActor(actor.second);
    }
    void PoseRenderer::removeAllActors()
    {
        LOCK
        auto colors =
            vtkSmartPointer<vtkNamedColors>::New();
        window->GetRenderWindow()
                ->RemoveRenderer(renderer);
        renderer =
            vtkSmartPointer<vtkRenderer>::New();
        window->GetRenderWindow()->AddRenderer(renderer);
        renderer->SetBackground(colors->GetColor3d("AliceBlue").GetData());
         window->GetRenderWindow()->MakeCurrent();
        posesAndActors.clear();
    }

    template vtkSmartPointer<vtkActor>
    PoseRenderer::displayPoints(const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> &, const char *, double);
}
