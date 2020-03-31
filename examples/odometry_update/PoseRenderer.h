//
// Created by tomlucas on 23.03.20.
//

#ifndef CYLINDEREXAMPLE_POSERENDERER_H
#define CYLINDEREXAMPLE_POSERENDERER_H
#include <vtkActor.h>
#include <vtkCommand.h>
#include <vtkNamedColors.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkCubeSource.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <vtkAnimationScene.h>
#include <vtkAnimationCue.h>

#include <list>
#include <ADEKF/types/SO3.h>
#include <memory>
#include <chrono>
#include <thread>

namespace adekf::viz {

    class PoseReader{
    public:
        virtual Eigen::Vector3d getPosition() =0;
        virtual adekf::SO3d getOrientation()=0;
        virtual ~PoseReader(){};
    };



    template<class EstimatorType>
    class GenericPoseReader : public PoseReader{
        EstimatorType * estimator;
    public:
        GenericPoseReader(EstimatorType * estimator):estimator(estimator){

        }
        virtual ~GenericPoseReader(){};

        Vector3d getPosition() override {
            return estimator->mu.position;
        }

        SO3d getOrientation() override {
            return estimator->mu.orientation;
        }
    };
    template<class EstimatorType>
    class GenericPositionReader : public PoseReader{
        EstimatorType * estimator;
    public:
        GenericPositionReader(EstimatorType * estimator):estimator(estimator){

        }
        virtual ~GenericPositionReader(){}


        Vector3d getPosition() override {
            return estimator->mu.position;
        }

        SO3d getOrientation() override {
            return SO3d::Identity();
        }
    };

    class PoseRenderer {
        //the windowInteractor that contains all graphic objects
        inline static vtkSmartPointer<vtkRenderWindowInteractor> window;
        inline static vtkSmartPointer<vtkRenderer> renderer;
        inline static std::list<std::pair<PoseReader *,vtkSmartPointer<vtkActor>  >> posesAndActors;


    public:

        static bool isDone(){
            return window->GetDone();
        }
        /**
         * Initializes the graphical interface
         */
        static void initGui() {
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
            renderWindow->SetSize(800,800);
            renderer->ResetCamera();
            renderer->GetActiveCamera()->SetPosition(0,0,80);
            renderer->GetActiveCamera()->SetFocalPoint(20,0,0);

            // Add the actor to the scene
            renderer->SetBackground(colors->GetColor3d("AliceBlue").GetData());
           renderWindow->MakeCurrent();
            // Render and interact

            vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
                    vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview
            window->SetInteractorStyle( style );
            window->Initialize();
            vtkSmartPointer<vtkCallbackCommand> timerCallback =
                    vtkSmartPointer<vtkCallbackCommand>::New();
            timerCallback->SetCallback ( [](auto a, auto eventId, auto c, auto d){
                if(eventId==vtkCommand::TimerEvent)
                    std::cout << "Got timer Event" << std::endl;
                if(eventId==vtkCommand::ExitEvent){
                    std::cout << "Got Exit Event" << std::endl;
                }
                std::cout << "Catched an Event" << std::endl;
                window->Render();

            } );

            window->AddObserver ( vtkCommand::MouseMoveEvent, timerCallback );
            renderWindow->Render();
        }

        static void updateWindow(){
            for(auto poseActor : posesAndActors){
                auto actor= poseActor.second;
                auto pose=poseActor.first;
                auto position =pose->getPosition();
                actor->SetPosition(position.x(),position.y(),position.z());
                auto q=pose->getOrientation();
                actor->SetOrientation(0,0,0);
                Eigen::AngleAxis<decltype(q)::Scalar> aa(q);
                auto axis=aa.axis();
                actor->RotateWXYZ(aa.angle()*180./M_PI,axis.x(),axis.y(),axis.z());
            }
            window->ProcessEvents();
            window->Render();
        }



        template<class EstimatorType>
        static void displayPosition(EstimatorType * estimator, const char * color){
            displayPoseGeneric<EstimatorType,GenericPositionReader<EstimatorType> >(estimator,color);
        }
        template<class EstimatorType>
        static void displayPose(EstimatorType * estimator, const char * color){
            displayPoseGeneric<EstimatorType,GenericPoseReader<EstimatorType> >(estimator,color);
        }

        template<class EstimatorType, class ReaderType>
        static void displayPoseGeneric(EstimatorType *estimator, const char * color) {
            auto colors =
                    vtkSmartPointer<vtkNamedColors>::New();
            // Create a sphere
            auto cubeSource =
                    vtkSmartPointer<vtkCubeSource>::New();
            cubeSource->SetCenter(0.0, 0.0, 0.0);
            cubeSource->SetBounds(-1,1,-0.5,0.5,-0.2,0.2);
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
            posesAndActors.push_back(std::make_pair<PoseReader *, vtkSmartPointer<vtkActor> >(new ReaderType(estimator),std::move(actor)));
        }
        static void disposeWindow(){
            for(auto poseActor : posesAndActors){
                delete poseActor.first;
            }
        }
    };

}


#endif //CYLINDEREXAMPLE_POSERENDERER_H
