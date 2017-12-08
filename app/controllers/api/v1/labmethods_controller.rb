class Api::V1::LabmethodsController < ApiController

def index
  @labmethods = Labmethod.all
  render json: @labmethods
end

def show
  # binding.pry
  selected_method = params[:id]
  filtered_method =Methodtype.find_by(method_type: selected_method)
  @labmethods = filtered_method.labmethods
  # binding.pry
  render json: @labmethods
end

end
